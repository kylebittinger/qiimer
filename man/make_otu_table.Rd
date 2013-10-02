\name{make_otu_table}
\alias{make_otu_table}
\title{Create an OTU table.}
\usage{
  make_otu_table(otus, sample_ids = NULL)
}
\arguments{
  \item{otus}{A list of QIIME-format sequence identifiers
  for each OTU.}

  \item{sample_ids}{An optional vector of sample IDs to
  include in the result.}
}
\value{
  A matrix of OTU counts by sample.
}
\description{
  Create an OTU table.
}

