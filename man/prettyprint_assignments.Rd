\name{prettyprint_assignments}
\alias{prettyprint_assignments}
\title{Reformat taxonomic assignments for presentation.}
\usage{
  prettyprint_assignments(assignments, sep = ";",
  prefix_rank = 3, suffix_rank = 6)
}
\arguments{
  \item{assignments}{A character vector of taxonomic
  assignments.}

  \item{sep}{The character separating taxa in each
  assignment.}

  \item{prefix_rank}{The rank of taxonomy to use as the
  first word in the prettyprinted assignment.}

  \item{suffix_rank}{The rank of taxonomy to use as the
  second word in the prettyprinted assignment.}
}
\value{
  A character vector of reformatted assignments.
}
\description{
  Reformat taxonomic assignments for presentation.
}

