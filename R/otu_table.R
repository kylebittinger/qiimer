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
prettyprint_assignments <- function(assignments, sep="; ", prefix_rank=2, suffix_rank=100) {
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

#' Standard taxonomic ranks
#'
#' @export
taxonomic_ranks <- c(
  "Domain", "Kingdom", "Phylum", "Class", "Order", 
  "Family", "Genus", "Species", "Strain")

#' Parse taxonomic assignment strings
#'
#' @param assignments Character vector of taxonomic assignments.
#' @param split Pattern on which to split taxa in assignment strings.
#' @param ranks Character vector of taxonomic ranks, used as row names in the
#'   resultant data frame.
#' @param ... Additional parameters are passed to the \code{strsplit} function.
#' @return A data frame of taxonomic assignments.
#' @export
parse_assignments <- function(
  assignments, split="; ", ranks=taxonomic_ranks[2:8], ...) {
  a <- strsplit(as.character(assignments), split, ...)
  max_ranks <- max(sapply(a, length))
  a <- lapply(a, function (x) {
    fill_length <- max_ranks - length(x)
    c(x, rep(NA, fill_length))
  })
  a <- as.data.frame(do.call(rbind, a))
  colnames(a) <- ranks[1:ncol(a)]
  a
}

#' Create a heatmap of OTU counts.
#'
#' @param otu_counts A matrix of OTU counts, one row per OTU and one column 
#'   per sample.
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts.
#' @param threshold Minimum number of OTU counts necessary for an assignment to
#'   be included in the heatmap.  Assignment groups are filtered prior to
#'   calculating the proportions, so these groups are effectively removed from
#'   the analysis.
#' @param color Vector of colors to use in the heatmap.
#' @param ... Additional arguments are passed to the pheatmap function.
#' @return A heatmap plot of the proportions of assignments in each sample
#' @export
otu_heatmap <- function(otu_counts, assignments, threshold=0,
  color=saturated_rainbow(max(colSums(otu_counts))),
  breaks=seq(0, 1, length.out=length(color) + 1), ...) {
  # The rowsum() function does not play well with factors, convert to character.
  assignments <- as.character(assignments)
  # NA values in assignments will work fine but produce a warning from rowsum().
  # However, NA values work well for representing unassigned taxa.  We therefore
  # suppress warnings from this function.
  assignment_counts <- suppressWarnings(rowsum(otu_counts, assignments))
  rows_to_keep <- apply(assignment_counts, 1, sum) >= threshold
  assignment_fracs <- apply(assignment_counts, 2, function (x) {x / sum(x)}) 
  assignment_fracs <- assignment_fracs[rows_to_keep,]
  pheatmap(as.matrix(assignment_fracs), color=color, breaks=breaks, ...)
}

#' Saturated rainbow palette
#'
#' @param n Length of the palette
#' @param saturation_limit The fraction of the total palette length over which
#'   the rainbow extends.  Above this limit, the color will remain the same.
#' @return A vector of colors.
#' @export
saturated_rainbow <- function (n, saturation_limit=0.4) {
  saturated_len <- floor(n * (1 - saturation_limit))
  unsaturated_len <- n - saturated_len - 2
  
  rainbow_colors <- rev(rainbow(unsaturated_len, start=0, end=0.6))
  
  first_color <- rainbow_colors[1]
  presence_colors <- colorRampPalette(c("#FFFFFFFF", first_color))(2)
  
  last_color <- rainbow_colors[length(rainbow_colors)]
  saturated_colors <- rep(last_color, saturated_len)
  
  c(presence_colors, rainbow_colors, saturated_colors)
}

#' Create a barplot of OTU assignments.
#'
#' @param otu_counts Matrix of OTU counts, one row per OTU and one column 
#'   per sample
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts.
#' @param max_categories  Maximum number of different assignments to be 
#'   displayed in the plot.  The top assignments are determined by total count, 
#'   and others are labeled "Other."  The default, 8, is the maximum supported 
#'   by palettes in color brewer.
#' @return A bar plot of OTU assignments.
#' @export
otu_barplot <- function(otu_counts, assignments, max_categories=8) {
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
    opts(
      axis.text.x=(theme_text(angle=90, hjust=1)), 
      panel.grid.major=theme_blank(),
      panel.grid_minor=theme_blank(),
      panel.border=theme_blank())
}
