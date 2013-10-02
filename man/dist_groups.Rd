\name{dist_groups}
\alias{dist_groups}
\title{Create a data frame of distances between groups of items.}
\usage{
  dist_groups(d, g)
}
\arguments{
  \item{d}{A distance matrix object of class `"dist"`.}

  \item{g}{A factor representing the groups of objects in
  `d`.}
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
\examples{
data(relmbeta_dist)
data(relmbeta)
dist_groups(relmbeta_dist, relmbeta$Diet)
}

