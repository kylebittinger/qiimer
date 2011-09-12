\name{parse_rarefaction}
\alias{parse_rarefaction}
\title{Read a collated alpha diversity table from QIIME.}
\usage{parse_rarefaction(filepath)}
\description{Read a collated alpha diversity table from QIIME.}
\value{A data frame representing the table.  The table is returned in
melted form so sample IDs appear in rows rather than in columns.  The
column of diversity values is named "diversity".}
\arguments{\item{filepath}{Path to alpha diversity table file.}}
