\name{parse_assignments}
\alias{parse_assignments}
\title{Parse taxonomic assignment strings}
\usage{
  parse_assignments(assignments, split = "; ",
    ranks = taxonomic_ranks[2:8], ...)
}
\arguments{
  \item{assignments}{Character vector of taxonomic
  assignments.}

  \item{split}{Pattern on which to split taxa in assignment
  strings.}

  \item{ranks}{Character vector of taxonomic ranks, used as
  row names in the resultant data frame.}

  \item{...}{Additional parameters are passed to the
  \code{strsplit} function.}
}
\value{
  A data frame of taxonomic assignments.
}
\description{
  Parse taxonomic assignment strings
}

