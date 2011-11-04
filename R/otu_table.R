#' Parse a QIIME OTU table file.
#'
#' @param filepath Path to OTU table file.
#' @param commented TRUE if the header line is preceeded by an additional
#'   comment line, otherwise FALSE.  This is usually the case for OTU
#'   tables generated with QIIME, so we default to TRUE.
#' @return A list with four attributes: sample_ids, otu_ids, counts, and 
#'   metadata, a data structure similar to that returned by the python 
#'   function qiime.parse.parse_otu_table.  The sample_ids, otu_ids, and
#'   metadata attributes are character vectors.  The counts attribute is an
#'   integer matrix with one column per sample_id and one row per otu_id.
#' @export
parse_otu_table <- function(filepath, commented=TRUE) {
  f <- file(filepath, "rt")
  header_line <- readLines(f, n=1)
  if (commented) {
    header_line <- readLines(f, n=1)
  }
  col_names <- unlist(strsplit(header_line, "\t"))

  col_classes <- rep("numeric", times=length(col_names))
  col_classes[c(1, length(col_classes))] <- "character"

  full_otu_table <- read.table(
    f, col.names=col_names, colClasses=col_classes, sep="\t", 
    quote="", as.is=TRUE, header=FALSE)
  close(f)

  data_cols <- 2:(length(col_names) - 1)

  sample_ids <- col_names[data_cols]
  otu_ids <- as.character(full_otu_table[,1])

  counts <- full_otu_table[,data_cols]
  rownames(counts) <- otu_ids

  metadata <- full_otu_table[,length(col_names)]
  names(metadata) <- otu_ids

  list(
    sample_ids=sample_ids, otu_ids=otu_ids, counts=counts, metadata=metadata)
}

#' Reformat taxonomic assignments for presentation.
#'
#' @param assignments A character vector of taxonomic assignments.
#' @param sep The character separating taxa in each assignment.
#' @param prefix_rank The rank of taxonomy to use as the first word in the
#'   prettyprinted assignment.
#' @param suffix_rank The rank of taxonomy to use as the second word in the
#'   prettyprinted assignment.
#' @return A character vector of reformatted assignments.
#' @export
prettyprint_assignments <- function(assignments, sep=";", prefix_rank=3, suffix_rank=6) {
  # make list of taxa
  alltaxa <- strsplit(assignments, sep)
  # Phylum and final classification
  superraw.labels <- lapply(alltaxa, function(x) {
    nx <- length(x)
    if (nx >= suffix_rank) x[c(prefix_rank, suffix_rank)] else
      if (nx > prefix_rank) x[c(prefix_rank, nx)] else x
    })
  # Paste these two items into a single string
  raw.labels <- lapply(superraw.labels, paste, collapse=" ")
  # Convert resulting nested list into a vector
  unlist(raw.labels)
}

#' Create a heatmap of OTU counts.
#'
#' @param otu_counts A matrix of OTU counts, one row per OTU and one column 
#'   per sample
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts
#' @param threshold Minimum number of OTU counts necessary for an assignment to
#'   be included in the heatmap.  Assignment groups are filtered prior to
#'   calculating the proportions, so these groups are effectively removed from
#'   the analysis.
#' @param color Vector of colors to use in the heatmap.  The default value is
#'   based on a 7-level GnBu sequential color scale from
#'   \url{http://colorbrewer2.org/}.  This differs from the default value in
#'   the pheatmap package, which can be found in the pheatmap documentation.
#' @param breaks Optional sequence of numbers used to map heatmap values to
#'   colors.  Must cover the range from 0 to 1 and be one element longer than
#'   the color vector. The default breaks are adapted to show relative
#'   abundance of major groups, while maintaining some contrast at the
#'   presence/absence level.  If value is NA (default for pheatmap), then the
#'   breaks are calculated automatically.
#' @param ... Additional arguments are passed to the pheatmap function.
#' @return A heatmap plot of the proportions of assignments in each sample
#' @export
otu_heatmap <- function(otu_counts, assignments, threshold=0,
  color=brewer.pal(7, "GnBu"), 
  breaks=c(0, 1e-05, 0.001, 0.01, 0.1, 0.2, 0.3, 1), ...) {
  assignment_counts <- rowsum(otu_counts, assignments)
  rows_to_keep <- apply(assignment_counts, 1, sum) >= threshold
  assignment_fracs <- apply(assignment_counts, 2, function (x) {x / sum(x)}) 
  assignment_fracs <- assignment_fracs[rows_to_keep,]
  pheatmap(as.matrix(assignment_fracs), breaks=breaks, color=color, ...)
}

#' Create a barplot of OTU assignments.
#'
#' @param otu_counts Matrix of OTU counts, one row per OTU and one column 
#'   per sample
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts
#' @return A bar plot of OTU assignments.
#' @export
otu_barplot <- function(otu_counts, assignments) {
  max_categories <- 8 # Max. supported by color brewer
  assignment_counts <- rowsum(otu_counts, assignments)
  assignment_labels <- rownames(assignment_counts)
  if (length(assignments) > max_categories) {
    assignment_order <- order(apply(assignment_counts, 1, sum), decreasing=TRUE)
    top_assignments <- assignment_order[1:max_categories]
    assignment_labels[-top_assignments] <- "Other"
    assignment_counts <- rowsum(assignment_counts, assignment_labels)
    assignment_labels <- rownames(assignment_counts)
  }
  melted_counts <- melt(
    data.frame(Assignment=assignment_labels, assignment_counts),
    variable_name="sample_id")
  ggplot(melted_counts, aes(x=sample_id, y=value, fill=Assignment)) + 
    geom_bar(position="fill") + 
    theme_bw() + 
    scale_x_discrete(name="", expand=c(0, 0)) + 
    scale_y_continuous(name="Percent composition", formatter="percent", expand=c(0, 0)) +
    scale_fill_brewer(palette="Set1") + 
    opts(
      axis.text.x=(theme_text(angle=90, hjust=1)), 
      panel.grid.major=theme_blank(),
      panel.grid_minor=theme_blank(),
      panel.border=theme_blank())
}
