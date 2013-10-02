\name{biom_raw_data}
\alias{biom_raw_data}
\title{Extract biom data in raw form.}
\usage{
  biom_raw_data(b)
}
\arguments{
  \item{b}{An object of class \code{"biom"}, typically
  created by the \code{read_biom} function in the
  \code{biom} library.}
}
\value{
  For sparse biom objects, returns a data frame.  For dense
  biom objects, returns a matrix.
}
\description{
  Extract biom data in raw form.
}
\examples{
\dontrun{
data(relmbeta_biom)
head(biom_raw_data(relmbeta_biom))
}
}

