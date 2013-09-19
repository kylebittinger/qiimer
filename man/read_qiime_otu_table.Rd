\name{read_qiime_otu_table}
\alias{read_qiime_otu_table}
\title{Parse a QIIME OTU table file in "calssic" format.}
\usage{
  read_qiime_otu_table(filepath, commented = TRUE)
}
\arguments{
  \item{filepath}{Path to OTU table file.}

  \item{commented}{TRUE if the header line is preceeded by
  an additional comment line, otherwise FALSE.  This is
  usually the case for OTU tables generated with QIIME, so
  we default to TRUE.}
}
\value{
  A list with four attributes: sample_ids, otu_ids, counts,
  and metadata, a data structure similar to that returned
  by the python function `qiime.parse.parse_otu_table`.
  The sample_ids, otu_ids, and metadata attributes are
  character vectors.  The counts attribute is an integer
  matrix with one column per sample_id and one row per
  otu_id.
}
\description{
  Parse a QIIME OTU table file in "calssic" format.
}

