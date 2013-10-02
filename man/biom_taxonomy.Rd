\name{biom_taxonomy}
\alias{biom_taxonomy}
\title{Extract taxonomy info from a biom object.}
\usage{
  biom_taxonomy(b, attr = "taxonomy")
}
\arguments{
  \item{b}{An object of class \code{"biom"}, typically
  created by the \code{read_biom} function in the `biom`
  library.}

  \item{attr}{The metadata attribute under which the
  taxonomy information can be found for each row item in
  the biom file.}
}
\value{
  A list of character vectors, one per row.
}
\description{
  Extract taxonomy info from a biom object.
}
\examples{
\dontrun{
data(relmbeta_biom)
head(biom_taxonomy(relmbeta_biom))
}
}

