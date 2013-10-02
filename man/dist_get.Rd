\name{dist_get}
\alias{dist_get}
\title{Retrieve distances from a `"dist"` object.}
\usage{
  dist_get(d, idx1, idx2)
}
\arguments{
  \item{d}{A distance matrix object of class `"dist"`.}

  \item{idx1,idx2}{Indicies specifying the distances to
  extract.}
}
\value{
  A vector of distances.
}
\description{
  Retrieve distances from a `"dist"` object.
}
\examples{
data(relmbeta_dist)
dist_get(relmbeta_dist, "A1", "A2")
dist_get(relmbeta_dist, "A1", c("A2", "A3", "A4", "A5"))
a <- paste("A", 1:5, sep="")
b <- paste("B", 1:5, sep="")
dist_get(relmbeta_dist, a, b)
}

