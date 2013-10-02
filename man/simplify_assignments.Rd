\name{simplify_assignments}
\alias{simplify_assignments}
\title{Reformat taxonomic assignments for presentation.}
\usage{
  simplify_assignments(assignments_df, rank1 = "Phylum",
    rank2 = "Genus")
}
\arguments{
  \item{assignments_df}{A data frame of taxonomic
  assignments.}

  \item{rank1}{The rank of taxonomy to use as the first
  word in the label.}

  \item{rank2}{The rank of taxonomy to use as the second
  word in the label.}
}
\value{
  A character vector of reformatted assignment labels.
}
\description{
  Reformat taxonomic assignments for presentation.
}
\examples{
data(relmbeta_assignments)
a <- split_assignments(relmbeta_assignments)
simplify_assignments(a)
}
\seealso{
  \code{\link{split_assignments}}
}

