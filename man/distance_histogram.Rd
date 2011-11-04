\name{distance_histogram}
\alias{distance_histogram}
\title{Plot a histogram of within-group and between-group distances.}
\usage{
  distance_histogram(distmat, sample_ids1, sample_ids2,
  breaks = 12, freq = TRUE, col1 = "blue", col2 = rgb(0, 1,
  0, 0.5), xlim = c(0, 1), main = "Distance Histogram",
  xlab = "Distance", legend_labels = c("Within-group",
  "Between-group"), ...)
}
\arguments{
  \item{distmat}{A square matrix of sample-sample
  distances, with sample IDs as row and column labels.}

  \item{sample_ids1}{Sample IDs in Group 1}

  \item{sample_ids2}{Sample IDs in Group 2}

  \item{breaks}{Number of bars in the histogram.}

  \item{freq}{If TRUE, plot absolute counts, if FALSE, plot
  relative density.}

  \item{col1}{Color of Group 1 in plot.}

  \item{col2}{Color of Group 2 in plot.}

  \item{xlim}{Limits of x-axis for histogram plot.}

  \item{main}{Title of plot.}

  \item{xlab}{Label of x-axis.}

  \item{legend_labels}{Character vector of length 2, giving
  labels for the legend.}

  \item{...}{passed directly to the plot() function.}
}
\description{
  Plot a histogram of within-group and between-group
  distances.
}

