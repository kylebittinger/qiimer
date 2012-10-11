\name{paired_distances}
\alias{paired_distances}
\title{Create a data frame of distances between paired samples}
\usage{
  paired_distances(distmat, sample_ids1, sample_ids2)
}
\arguments{
  \item{distmat}{A square matrix of sample-sample
  distances, with sample IDs as row and column labels.}

  \item{sample_ids1}{Sample IDs in Group 1.}

  \item{sample_ids2}{Sample IDs in Group 2.  The sample IDs
  in sample_ids2 will be paired with the sample IDs in
  sample_ids1 using cbind().}
}
\value{
  A data frame with 4 columns. "SampleA" and "SampleB"
  identify the samples compared. "Category" is a factor
  column indicating the type of comparison. Its values are:
  "Between Pairs" for upaired samples and "Within Pairs"
  for paired samples. "Distance" contains the applicable
  entry from the distance matrix.
}
\description{
  Create a data frame of distances between paired samples
}

