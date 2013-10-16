\name{biom_raw_data}
\alias{biom_raw_data}
\title{Extract raw data from a BIOM object.}
\usage{
  biom_raw_data(b)
}
\arguments{
  \item{b}{A BIOM object.}
}
\value{
  For sparse biom objects, returns a 3-column data frame of
  row names, column names, and the data value.  The first
  column is named using the first word in the BIOM object's
  \code{type} attribute (e.g. "OTU" for OTU tables).  The
  second and third columns are named "SampleID" and
  "value", respectively.  For dense biom objects, returns a
  matrix.
}
\description{
  This function extracts the raw data from a BIOM object,
  using the correct row and column names in the result.
}
\details{
  The BIOM object can be any list-like representation of
  the JSON source code in a BIOM-format file produced by
  QIIME.  There are several options for creating BIOM
  objects from QIIME output files.  The official library
  for BIOM files, \code{biom}, can create compatible
  objects via the \code{read_biom} function.  Alternately,
  the \code{fromJSON} function from either \code{RJSONIO}
  or \code{rjson} may be used.
}
\examples{
data(relmbeta_biom)
head(biom_raw_data(relmbeta_biom))
}

