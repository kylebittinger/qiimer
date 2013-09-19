\name{dist_groups}
\alias{dist_groups}
\title{Create a data frame of distances between groups of items.}
\usage{
  dist_groups(d, g)
}
\arguments{
  \item{d}{An object of class `dist`, representing a
  pairwise distance matrix.  A square matrix can also be
  passed; it will be converted to a `dist` object using
  `as.dist`.}

  \item{g}{A factor representing the groups of objects in
  the matrix.}
}
\value{
  A data frame with 6 columns. "Item1" and "Item2" identify
  the items compared, using the label if available.
  Likewise, "Group1" and "Group2" identify the groups of
  the items. "Label" is a factor giving a convenient label
  for the type of comparison. Finally, "Distance" contains
  the distance of interest.
}
\description{
  Create a data frame of distances between groups of items.
}

