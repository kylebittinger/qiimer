\name{otu_barplot}
\alias{otu_barplot}
\title{Create a barplot of OTU assignments.}
\usage{
  otu_barplot(otu_counts, assignments)
}
\arguments{
  \item{otu_counts}{Matrix of OTU counts, one row per OTU
  and one column per sample}

  \item{assignments}{A character vector of OTU assignments.
  Length should match number of rows in otu_counts}
}
\value{
  A bar plot of OTU assignments.
}
\description{
  Create a barplot of OTU assignments.
}

