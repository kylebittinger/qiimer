#' Read a distance matrix file and return a matrix.
#'
#' @param filepath Path to distance matrix file.
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
