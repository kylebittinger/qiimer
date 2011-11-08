\name{grouped_distances}
\alias{grouped_distances}
\title{Create a data frame of distances between groups of samples.}
\usage{
  grouped_distances(dm, sample_ids1, sample_ids2)
}
\arguments{
  \item{distmat}{A square matrix of sample-sample
  distances, with sample IDs as row and column labels.}

  \item{sample_ids1}{Sample IDs in Group 1}

  \item{sample_ids2}{Sample IDs in Group 2}
}
\value{
  A data frame with 4 columns. "SampleA" and "SampleB"
  identify the samples compared. "Category" is a factor
  column indicating the type of comparison. Its values are:
  "Between" for comparisons between sample_ids1 and
  sample_ids2, "Within1" for comparisons within
  sample_ids1, and "Within2" for comparisons within
  sample_ids2. "Distance" contains the applicable entry
  from the distance matrix.
}
\description{
  Create a data frame of distances between groups of
  samples.
}

