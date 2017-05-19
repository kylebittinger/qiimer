#' Read a QIIME distance matrix file.
#'
#' @param filepath Path to QIIME distance matrix file.
#' @return A distance matrix.
#' @export
read_qiime_distmat <- function(filepath) {
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
  as.dist(distmat)
}

#' Retrieve distances from a `"dist"` object.
#' 
#' DEPRECATED: This function has been moved to the usedist package.
#'
#' @param d A distance matrix object of class `"dist"`.
#' @param idx1,idx2 Indicies specifying the distances to extract.
#' @return A vector of distances.
#' @export
#' @examples
#' data(relmbeta_dist)
#' dist_get(relmbeta_dist, "A1", "A2")
#' dist_get(relmbeta_dist, "A1", c("A2", "A3", "A4", "A5"))
#' dist_get(relmbeta_dist, c("A1", "A2", "A3"), c("B1", "B2", "B3"))
dist_get <- function (d, idx1, idx2) {
  warning(
    "Deprecated: This function has been moved to the usedist package.",
    call. = FALSE)
  usedist::dist_get(d, idx1, idx2)
}

#' Extract parts of a `"dist"` object.
#'
#' DEPRECATED: This function has been moved to the usedist package.
#'
#' @param d A distance matrix object of class `"dist"`.
#' @param idx Indices specifying the subset of distances to extract.
#' @return A distance matrix.
#' @export
#' @examples
#' data(relmbeta_dist)
#' dist_subset(relmbeta_dist, c("A1", "A2", "A3", "A4", "A5"))
dist_subset <- function (d, idx) {
  warning(
    "Deprecated: This function has been moved to the usedist package.",
    call. = FALSE)
  usedist::dist_subset(d, idx)
}

#' Create a data frame of distances between groups of items.
#'
#' DEPRECATED: This function has been moved to the usedist package.
#'
#' @param d A distance matrix object of class `"dist"`.
#' @param g A factor representing the groups of objects in `d`.
#' @return A data frame with 6 columns. "Item1" and "Item2" identify the
#'   items compared, using the label if available. Likewise, "Group1" and 
#'   "Group2" identify the groups of the items. "Label" is a factor giving a
#'   convenient label for the type of comparison. Finally, "Distance" contains
#'   the distance of interest.
#' @export
#' @examples
#' data(relmbeta_dist)
#' data(relmbeta)
#' head(dist_groups(relmbeta_dist, relmbeta$Diet))
dist_groups <- function(d, g) {
  warning(
    "Deprecated: This function has been moved to the usedist package.",
    call. = FALSE)
  usedist::dist_groups(d, g)
}
