#' Read a distance matrix file.
#'
#' @param filepath Path to distance matrix file.
#'
#' @return A matrix of sample-to-sample distances.
ReadDistanceMatrix <- function(filepath) {
  distmat.file <- file(filepath, 'rt')
  header <- readLines(distmat.file, n=1)
  cols <- unlist(strsplit(header, "\t"))

  # Provide column classes to speed up conversion.
  col.classes <- rep("numeric", times=length(cols))
  # From the read.table docs: Note that colClasses is specified per
  # column (not per variable) and so includes the column of row names.
  col.classes[1] <- "character"
  
  as.matrix(read.table(distmat.file,
                       col.names=cols,
                       row.names=1,
                       colClasses=col.classes,
                       sep="\t",
                       quote=""))
}

#' Compare within-group to between-group distances.
GroupDistances <- function(distmat, sample.ids1, sample.ids2) {
  between <- distmat[sample.ids1, sample.ids2]
  dim(between) <- NULL # Flatten matrix
  within1 <- as.vector(as.dist(distmat[sample.ids1, sample.ids1]))
  within2 <- as.vector(as.dist(distmat[sample.ids2, sample.ids2]))
  list(between=between, within=c(within1, within2), within1=within1, within2=within2)
}

#' Plot a histogram of within-group and between-group distances.
DistanceHistogram <- function(distmat, sample.ids1, sample.ids2, breaks=12,
                              freq=TRUE, col1="blue", col2=rgb(0,1,0,0.5),
                              xlim=c(0,1),
                              main="Distance Histogram", xlab="Distance",
                              ...) {
  group.distances <- GroupDistances(distmat, sample.ids1, sample.ids2)
  #col2.a <- col2rgb(col2)
  #col2.b <- rgb(col2.a[1,1], col2.a[2,1], col2.a[3,1], alpha, maxColorValue=255)
  h1 <- hist(group.distances$within, breaks=breaks, plot=FALSE)
  h2 <- hist(group.distances$between, breaks=breaks, plot=FALSE)
  ymax <- if (freq) max(c(h1$counts, h2$counts)) else max(c(h1$density, h2$density))
  plot(h1, col=col1, freq=freq, ylim=c(0, ymax), xlim=xlim, main=main, xlab=xlab, ...)
  plot(h2, col=col2, freq=freq, add=T)
  legend("topleft", c("Within-group", "Between-group"), fill=c(col1, col2))
}
