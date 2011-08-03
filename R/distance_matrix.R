#' Read a distance matrix file.
#'
#' @param filepath Path to distance matrix file.
#'
#' @return A matrix of sample-to-sample distances.
parse_distmat <- function(filepath) {
  distmat_file <- file(filepath, 'rt')
  header_line <- readLines(distmat_file, n=1)
  column_names <- unlist(strsplit(header_line, "\t"))

  # Provide column classes to speed up conversion.
  column_classes <- rep("numeric", times=length(column_names))
  # From the read.table docs: Note that colClasses is specified per
  # column (not per variable) and so includes the column of row names.
  column_classes[1] <- "character"
  
  distmat <- as.matrix(read.table(
    distmat_file, col.names=column_names, row.names=1,
    colClasses=column_classes, sep="\t", quote=""))

  close(distmat_file)
  distmat
}

#' Compare within-group to between-group distances.
#'
#' @param distmat A square matrix of sample-sample distances, with sample IDs
#'   as row and column labels.
#'
#' @param sample_ids1 Sample IDs in Group 1
#'
#' @param sample_ids2 Sample IDs in Group 2
#'
#' @return A list with 4 attributes: between, within, within1, and within2.
#'   The attributes are assigned to vectors of distances between groups,
#'   within groups, within the first group, and within the second group,
#'   respectively.
group_distances <- function(distmat, sample_ids1, sample_ids2) {
  between <- distmat[sample_ids1, sample_ids2]
  dim(between) <- NULL # Flatten matrix
  within1 <- as.vector(as.dist(distmat[sample_ids1, sample_ids1]))
  within2 <- as.vector(as.dist(distmat[sample_ids2, sample_ids2]))
  list(between=between, within=c(within1, within2),
    within1=within1, within2=within2)
}

#' Plot a histogram of within-group and between-group distances.
#'
#' @param distmat A square matrix of sample-sample distances, with sample IDs
#'   as row and column labels.
#' @param sample_ids1 Sample IDs in Group 1
#' @param sample_ids2 Sample IDs in Group 2
#' @param breaks Number of bars in the histogram.
#' @param freq If TRUE, plot absolute counts, if FALSE, plot relative density.
#' @param col1 Color of Group 1 in plot.
#' @param col2 Color of Group 2 in plot.
#' @param xlim Limits of x-axis for histogram plot.
#' @param main Title of plot.
#' @param xlab Label of x-axis.
#' @param legend_labels Character vector of length 2, giving labels for the legend.
#' @param ... passed directly to the plot() function.
distance_histogram <- function(distmat, sample_ids1, sample_ids2, breaks=12,
  freq=TRUE, col1="blue", col2=rgb(0,1,0,0.5), xlim=c(0,1),
  main="Distance Histogram", xlab="Distance",
  legend_labels=c("Within-group", "Between-group"), ...) {
  group.distances <- GroupDistances(distmat, sample_ids1, sample_ids2)
  h1 <- hist(group.distances$within, breaks=breaks, plot=FALSE)
  h2 <- hist(group.distances$between, breaks=breaks, plot=FALSE)
  ymax <- if (freq) max(c(h1$counts, h2$counts)) else max(c(h1$density, h2$density))
  plot(h1, col=col1, freq=freq, ylim=c(0, ymax), xlim=xlim, main=main, xlab=xlab, ...)
  plot(h2, col=col2, freq=freq, add=T)
  legend("topleft", legend_labels, fill=c(col1, col2))
}
