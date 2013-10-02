\name{otu_heatmap}
\alias{otu_heatmap}
\title{Create a heatmap of OTU counts.}
\usage{
  otu_heatmap(otu_counts, assignments, threshold = 0,
    plot = T,
    color = saturated_rainbow(max(colSums(otu_counts)) + 1),
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
  Assignments are filtered after calculating the
  proportions, so the threshold setting does not affect the
  display of the remaining OTUs.}

  \item{plot}{If true, display a plot.  If false, just
  return the computed abundances.}

  \item{color}{Vector of colors to use in the heatmap.}

  \item{breaks}{Vector of color breaks, one element greater
  in length than `colors`.}

  \item{...}{Additional arguments are passed to the
  pheatmap function.}
}
\value{
  A heatmap plot of the proportions of assignments in each
  sample, and invisibly returns a matrix of the proportions
  in the plot.
}
\description{
  Create a heatmap of OTU counts.
}
\examples{
data(relmbeta_assignments)
data(relmbeta_counts)
a <- simplify_assignments(split_assignments(relmbeta_assignments))

\dontrun{
otu_heatmap(relmbeta_counts, a, threshold=10)
otu_heatmap(
  relmbeta_counts, a, threshold=10,
  cluster_rows=FALSE, cluster_cols=FALSE,
  cellwidth=12, cellheight=12)
}

heatmap_data <- otu_heatmap(relmbeta_counts, a, threshold=10, plot=FALSE)
head(heatmap_data)
}
\seealso{
  \code{\link{pheatmap}}, \code{\link{saturated_rainbow}}
}

