\name{group_distances}
\alias{group_distances}
\title{Compare within-group to between-group distances.}
\usage{
  group_distances(distmat, sample_ids1, sample_ids2)
}
\arguments{
  \item{distmat}{A square matrix of sample-sample
  distances, with sample IDs as row and column labels.}

  \item{sample_ids1}{Sample IDs in Group 1}

  \item{sample_ids2}{Sample IDs in Group 2}
}
\value{
  A list with 4 attributes: between, within, within1, and
  within2. The attributes are assigned to vectors of
  distances between groups, within groups, within the first
  group, and within the second group, respectively.
}
\description{
  Compare within-group to between-group distances.
}

