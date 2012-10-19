\name{otu_barplot}
\alias{otu_barplot}
\title{Create a barplot of OTU assignments.}
\usage{
  otu_barplot(otu_counts, assignments, max_categories = 8)
}
\arguments{
  \item{otu_counts}{Matrix of OTU counts, one row per OTU
  and one column per sample}

  \item{assignments}{A character vector of OTU assignments.
  Length should match number of rows in otu_counts.}

  \item{max_categories}{Maximum number of different
  assignments to be displayed in the plot.  The top
  assignments are determined by total count, and others are
  labeled "Other."  The default, 8, is the maximum
  supported by palettes in color brewer.}
}
\value{
  A bar plot of OTU assignments.
}
\description{
  Create a barplot of OTU assignments.
}

