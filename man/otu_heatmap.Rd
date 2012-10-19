\name{otu_heatmap}
\alias{otu_heatmap}
\title{Create a heatmap of OTU counts.}
\usage{
  otu_heatmap(otu_counts, assignments, threshold = 0,
    color = saturated_rainbow(max(colSums(otu_counts))),
    breaks = seq(0, 1, length.out = length(color) + 1),
    ...)
}
\arguments{
  \item{otu_counts}{A matrix of OTU counts, one row per OTU
  and one column per sample.}

  \item{assignments}{A character vector of OTU assignments.
  Length should match number of rows in otu_counts.}

  \item{threshold}{Minimum number of OTU counts necessary
  for an assignment to be included in the heatmap.
  Assignment groups are filtered prior to calculating the
  proportions, so these groups are effectively removed from
  the analysis.}

  \item{color}{Vector of colors to use in the heatmap.}

  \item{...}{Additional arguments are passed to the
  pheatmap function.}
}
\value{
  A heatmap plot of the proportions of assignments in each
  sample
}
\description{
  Create a heatmap of OTU counts.
}

