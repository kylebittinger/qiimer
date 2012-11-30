#' Read a distance matrix file.
#'
#' @param filepath Path to distance matrix file.
#' @return A matrix of sample-to-sample distances.
#' @export
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

#' Create a data frame of distances between groups of samples.
#'
#' @param distmat A square matrix of sample-sample distances, with sample IDs
#'   as row and column labels.
#' @param sample_ids1 Sample IDs in Group 1
#' @param sample_ids2 Sample IDs in Group 2
#' @return A data frame with 4 columns. "SampleA" and "SampleB" identify the
#'   samples compared. "Category" is a factor column indicating the type of
#'   comparison. Its values are: "Between" for comparisons between sample_ids1
#'   and sample_ids2, "Within1" for comparisons within sample_ids1, and
#'   "Within2" for comparisons within sample_ids2. "Distance" contains the
#'   applicable entry from the distance matrix.
#' @export
grouped_distances <- function (distmat, sample_ids1, sample_ids2) {
  df <- expand.grid(SampleA=sample_ids1, SampleB=sample_ids2)
  categories <- rep("Between", length(df$SampleA))
  
  if (length(sample_ids1) > 1) {
    combin1 <- combn(sample_ids1, 2, simplify=FALSE)
    within1 <- do.call(rbind, lapply(combin1, function (x) {
      data.frame(SampleA=x[1], SampleB=x[2])
    }))
    df <- rbind(within1, df)
    categories <- c(rep("Within1", length(within1$SampleA)), categories)
  }
  
  if (length(sample_ids2) > 1) {
    combin2 <- combn(sample_ids2, 2, simplify=FALSE)
    within2 <- do.call(rbind, lapply(combin2, function (x) {
      data.frame(SampleA=x[1], SampleB=x[2])
    }))
    df <- rbind(df, within2)
    categories <- c(categories, rep("Within2", length(within2$SampleA)))
  }
    
  df$Category <- factor(categories, levels=c("Within1", "Between", "Within2"))
  df$Distance <- apply(df, 1, function (x) { distmat[x[1], x[2]] })
  df
}

#' Create a data frame of distances between paired samples
#'
#' @param distmat A square matrix of sample-sample distances, with sample IDs
#'   as row and column labels.
#' @param sample_ids1 Sample IDs in Group 1.
#' @param sample_ids2 Sample IDs in Group 2.  The sample IDs in sample_ids2 will
#'   be paired with the sample IDs in sample_ids1 using cbind().
#' @return A data frame with 4 columns. "SampleA" and "SampleB" identify the
#'   samples compared. "Category" is a factor column indicating the type of
#'   comparison. Its values are: "Between Pairs" for upaired samples and 
#'   "Within Pairs" for paired samples. "Distance" contains the
#'   applicable entry from the distance matrix.
#' @export
paired_distances <- function(distmat, sample_ids1, sample_ids2) {
  df <- expand.grid(
    SampleA=sample_ids1,
    SampleB=sample_ids2)

  # Find row index of paired samples in the data frame
  # More efficient solutions rely on implementation of expand.grid:
  #   n <- length(sample_ids1)
  #   within_pair_idx <- seq(1, n * n, by=(n + 1))
  pairs <- cbind(sample_ids1, sample_ids2)
  within_pair_idx <- apply(pairs, 1, function(x) {
    which((df$SampleA == x[1]) & (df$SampleB == x[2]))
  })
  categories <- rep("Between Pairs", length(df$SampleA))
  categories[within_pair_idx] <- "Within Pairs"
  df$Category <- factor(categories, levels=c("Within Pairs", "Between Pairs"))

  df$Distance <- apply(df, 1, function (x) {
    distmat[x[1], x[2]]
  })
  df
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
#' @export
distance_histogram <- function(distmat, sample_ids1, sample_ids2, breaks=12,
  freq=TRUE, col1="blue", col2=rgb(0,1,0,0.5), xlim=c(0,1),
  main="Distance Histogram", xlab="Distance",
  legend_labels=c("Within-group", "Between-group"), ...) {
  group.distances <- grouped_distances(distmat, sample_ids1, sample_ids2)
  h1 <- hist(group.distances$within, breaks=breaks, plot=FALSE)
  h2 <- hist(group.distances$between, breaks=breaks, plot=FALSE)
  ymax <- if (freq) max(c(h1$counts, h2$counts)) else max(c(h1$density, h2$density))
  plot(h1, col=col1, freq=freq, ylim=c(0, ymax), xlim=xlim, main=main, xlab=xlab, ...)
  plot(h2, col=col2, freq=freq, add=T)
  legend("topleft", legend_labels, fill=c(col1, col2))
}
