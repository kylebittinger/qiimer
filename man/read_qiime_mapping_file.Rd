\name{read_qiime_mapping_file}
\alias{read_qiime_mapping_file}
\title{Read a QIIME sample mapping file.}
\usage{
  read_qiime_mapping_file(filepath)
}
\arguments{
  \item{filepath}{Path to sample mapping file.  The file
  must conform to the QIIME standards, detailed at
  \url{http://qiime.org/documentation/file_formats.html}.}
}
\value{
  A data frame of sample information.  Because the SampleID
  column is so often used to extract data from distance
  matrices and OTU tables, it is returned as a character
  vector.
}
\description{
  Read a QIIME sample mapping file.
}

