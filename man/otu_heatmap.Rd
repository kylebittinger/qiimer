\name{otu_heatmap}
\alias{otu_heatmap}
\title{Create a heatmap of OTU counts.}
\usage{
  otu_heatmap(otu_counts, assignments, threshold = 0, color
  = brewer.pal(7, "GnBu"), breaks = c(0, 1e-05, 0.001,
  0.01, 0.1, 0.2, 0.3, 1), ...)
}
\arguments{
  \item{otu_counts}{A matrix of OTU counts, one row per OTU
  and one column per sample}

  \item{assignments}{A character vector of OTU assignments.
  Length should match number of rows in otu_counts}

  \item{threshold}{Minimum number of OTU counts necessary
  for an assignment to be included in the heatmap.
  Assignment groups are filtered prior to calculating the
  proportions, so these groups are effectively removed from
  the analysis.}

  \item{color}{Vector of colors to use in the heatmap.  The
  default value is based on a 7-level GnBu sequential color
  scale from \url{http://colorbrewer2.org/}.  This differs
  from the default value in the pheatmap package, which can
  be found in the pheatmap documentation.}

  \item{breaks}{Optional sequence of numbers used to map
  heatmap values to colors.  Must cover the range from 0 to
  1 and be one element longer than the color vector. The
  default breaks are adapted to show relative abundance of
  major groups, while maintaining some contrast at the
  presence/absence level.  If value is NA (default for
  pheatmap), then the breaks are calculated automatically.}

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

