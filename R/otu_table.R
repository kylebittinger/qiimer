#' Parse a QIIME OTU table file in "calssic" format.
#'
#' @param filepath Path to OTU table file.
#' @param commented TRUE if the header line is preceeded by an additional
#'   comment line, otherwise FALSE.  This is usually the case for OTU
#'   tables generated with QIIME, so we default to TRUE.
#' @return A list with four attributes: sample_ids, otu_ids, counts, and 
#'   metadata, a data structure similar to that returned by the python 
#'   function `qiime.parse.parse_otu_table`.  The sample_ids, otu_ids, and
#'   metadata attributes are character vectors.  The counts attribute is an
#'   integer matrix with one column per sample_id and one row per otu_id.
#' @export
read_qiime_otu_table <- function(filepath, commented=TRUE) {
  f <- file(filepath, "rt")
  header_line <- readLines(f, n=1)
  if (commented) {
    header_line <- readLines(f, n=1)
  }
  col_names <- unlist(strsplit(header_line, "\t"))

  # Following the QIIME implementation, we decide if the table contains 
  # metadata by inspecting the final column for "Consensus Lineage", "OTU
  # Metadata", or "Taxonomy" (spaces and capitalization ignored)
  has_metadata <- grepl(
    "(Consensus ?Lineage)|(OTU ?Metadata)|(Taxonomy)",
    col_names[length(col_names)],
    ignore.case=T, perl=T)
  
  col_classes <- rep("numeric", times=length(col_names))
  
  if (has_metadata) {
    col_classes[c(1, length(col_classes))] <- "character"
  }
  
  full_otu_table <- read.table(
    f, col.names=col_names, colClasses=col_classes, sep="\t", 
    quote="", as.is=TRUE, header=FALSE)
  close(f)

  data_cols <- if (has_metadata) {
    2:(length(col_names) - 1) 
  } else {
    2:length(col_names)
  } 

  sample_ids <- col_names[data_cols]
  otu_ids <- as.character(full_otu_table[,1])

  counts <- as.matrix(full_otu_table[,data_cols])
  rownames(counts) <- otu_ids

  if (has_metadata) {
    metadata <- as.character(full_otu_table[,length(col_names)])
    names(metadata) <- otu_ids
  } else {
    metadata <- NULL
  }
    
  list(
    sample_ids=sample_ids, otu_ids=otu_ids, counts=counts, metadata=metadata)
}

#' Standard taxonomic ranks.
#'
#' @export
taxonomic_ranks <- c(
  "Domain", "Kingdom", "Phylum", "Class", "Order", 
  "Family", "Genus", "Species", "Strain")

#' Split taxonomic assignment strings
#'
#' @param assignments Character vector of taxonomic assignments.
#' @param ranks Character vector of taxonomic ranks, used as column names in the
#'   resultant data frame.
#' @param split Pattern on which to split taxa in assignment strings.
#' @param ... Additional parameters are passed to the \code{strsplit} function.
#' @return A data frame of taxonomic assignments.
#' @export
split_assignments <- function(assignments, ranks=taxonomic_ranks[2:8], 
  split="; ", ...) {
  a <- strsplit(as.character(assignments), split, ...)
  max_ranks <- max(sapply(a, length))
  a <- lapply(a, function (x) {
    fill_length <- max_ranks - length(x)
    c(x, rep(NA, fill_length))
  })
  a <- as.data.frame(do.call(rbind, a))
  colnames(a) <- ranks[1:ncol(a)]
  if (!is.null(names(assignments))) {
    rownames(a) <- names(assignments)
  }
  a
}

#' Reformat taxonomic assignments for presentation.
#'
#' @param assignments_df A data frame of taxonomic assignments.
#' @param rank1 The rank of taxonomy to use as the first word in the label.
#' @param rank2 The rank of taxonomy to use as the second word in the label.
#' @return A character vector of reformatted assignment labels.
#' @export
simplify_assignments <- function(assignments_df, rank1="Phylum", rank2="Genus") {
  if (is.character(rank1)) {
    rank1 <- match(rank1, colnames(assignments_df))
  }
  if (is.character(rank2)) {
    rank2 <- match(rank2, colnames(assignments_df))
  }
  apply(assignments_df, 1, function (x) {
    x <- na.omit(as.character(x))
    n <- length(x)
    if (n == 1)     return(x)
    if (n < rank1)  return(paste(x, collapse=" "))
    if (n == rank1) return(x[rank1])
    if (n < rank2)  return(paste(x[rank1], x[length(x)]))
    return(paste(x[rank1], x[rank2]))
  })
}

#' Create a heatmap of OTU counts.
#'
#' @param otu_counts A matrix of OTU counts, one row per OTU and one column 
#'   per sample.
#' @param assignments A character vector of OTU assignments.  Length should
#'   match number of rows in otu_counts.
#' @param threshold Minimum number of OTU counts necessary for an assignment to
#'   be included in the heatmap.  Assignments are filtered after calculating
#'   the proportions, so the threshold setting does not affect the display of
#'   the remaining OTUs.
#' @param plot If true, display a plot.  If false, just return the computed
#'   abundances.
#' @param color Vector of colors to use in the heatmap.
#' @param breaks Vector of color breaks, one element greater in length than
#'   `colors`.
#' @param ... Additional arguments are passed to the pheatmap function.
#' @return A heatmap plot of the proportions of assignments in each sample,
#'   and invisibly returns a matrix of the proportions in the plot.
#' @export
otu_heatmap <- function(otu_counts, assignments, threshold=0, plot=T,
  color=saturated_rainbow(max(colSums(otu_counts))),
  breaks=seq(0, 1, length.out=length(color) + 1), ...) {
  # rowsum() does not play well with factors
  assignments <- as.character(assignments)
  # NA values in assignments work fine, but produce an unnecessary warning
  assignment_counts <- suppressWarnings(rowsum(otu_counts, assignments))
  rows_to_keep <- apply(assignment_counts, 1, sum) >= threshold
  assignment_fracs <- apply(assignment_counts, 2, function (x) {x / sum(x)}) 
  assignment_fracs <- as.matrix(assignment_fracs[rows_to_keep,])
  if (plot) {
    pheatmap(assignment_fracs, color=color, breaks=breaks, ...)
  }
  # Return underlying data
  invisible(assignment_fracs)
}

#' Saturated rainbow palette.
#'
#' @param n Length of the palette
#' @param saturation_limit The fraction of the total palette length over which
#'   the rainbow extends.  Above this limit, the color will remain the same.
#' @return A vector of colors.
#' @export
saturated_rainbow <- function (n, saturation_limit=0.4) {
  saturated_len <- floor(n * (1 - saturation_limit))
  rainbow_colors <- rev(rainbow(n - saturated_len, start=0, end=0.6))
  last_color <- tail(rainbow_colors, n=1)
  saturated_colors <- rep(last_color, saturated_len)
  colors <- c(rainbow_colors, saturated_colors)
  colors[1] <- "#FFFFFFFF"
  colors
}
