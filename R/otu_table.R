#' Read a QIIME OTU table.
#'
#' @param filepath Path to OTU table file.
#'
#' @param commented TRUE if the header line is preceeded by an additional
#'   comment line, otherwise FALSE.  This is usually the case for OTU
#'   tables generated with QIIME, so we default to TRUE.
#'
#' @return A data frame representing the OTU table.  The data frame takes the
#'   form of an association table, where each row specifies an OTU, a SampleID,
#'   and a value for the count.  Assigned lineages are listed in a column named
#'   \code{Assignment}.  The format data frame is different from the
#'   text-formatted OTU table, which is a matrix of OTU's by Sample.
ReadOtuTable <- function(filepath, commented=TRUE) {
  otu.file <- file(filepath, 'rt')

  header <- readLines(otu.file, n=1)
  if (commented) {
    header <- readLines(otu.file, n=1)
  }

  cols <- unlist(strsplit(header, "\t"))
  cols[1] <- "OtuID"
  cols[length(cols)] <- "Assignment"
  
  otu.table <- read.table(otu.file,
                          col.names=cols,
                          sep="\t",
                          quote="",
                          as.is=TRUE,
                          row.names=NULL)
  close(otu.file)

  otu.table$OtuID <- as.factor(otu.table$OtuID)
  otu.table <- melt(otu.table,
                    id.vars=c('OtuID', 'Assignment'),
                    variable_name="SampleID",
                    )
  otu.table <- rename(otu.table, c(value="Counts"))
  otu.table
}

#' List OTU assignments.
#'
#' @param otu.table OTU table dataframe.
#'
#' @return A list of assignments per OTU.
SampleTotals <- function(otu.table) {
  tapply(otu.table$Assignment, otu.table$OtuID, function (x) { x[1] })
}

#' List total counts for each sample.
#'
#' @param otu.table OTU table dataframe.
#'
#' @return A list of total counts per sample.
SampleTotals <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$SampleID, sum)
}


#' List total counts for each OTU.
#'
#' @param otu.table OTU table dataframe.
#'
#' @return A list of total counts per OTU.
OtuTotals <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$OtuID, sum)
}


#' Tabulate OTU counts by sample.
#'
#' @param otu.table OTU table dataframe.
#'
#' @return A data frame with samples listed across the columns and OTU's
#'   appearing across the rows.  The form of this data frame resembles the
#'   original OTU table text file.
OtuCounts <- function(otu.table) {
  tapply(otu.table$Counts, list(otu.table$OtuID, otu.table$SampleID))
}


#' Tablulate counts per assigned lineage by sample.
#'
#' @param otu.table OTU table dataframe.
#'
#' @return A matrix with assigned lineages in the rows and samples in the
#'   columns.
LineageCounts <- function(otu.table) {
  tapply(otu.table$Counts, list(otu.table$Assignment, otu.table$SampleID), sum)
}

#' Attach sample metadata to OTU table data frame.
#'
#' @param otu.table OTU table dataframe.
#'
#' @param sample.mapping Sample mapping dataframe, containing a
#'   \code{SampleID} column.
AttachMetadata <- function(otu.table, sample.mapping) {
  merge(otu.table, sample.mapping, by=c('SampleID'))
}


#' Create a heatmap of taxonomic assignments.
#'
#' @param threshold Minimum number of OTU counts necessary for an assignment to
#'   be included in the heatmap.  Assignment groups are filtered prior to
#'   calculating the proportions, so these groups are effectively removed from
#'   the analysis.
#'
#' @param sample_order Optional vector of integers, specifying the order of the
#'   samples (columns) in the heatmap.  A suitable vector can be generated from
#'   a sample mapping table by applying the built-in \code{order} function to
#'   the column of interest.  By default, the samples are sorted in
#'   alphabetical order.
#'
#' @param color Vector of colors to use in the heatmap.  The default value is
#'   based on a 7-level GnBu sequential color scale from
#'   \url{http://colorbrewer2.org/}.  This differs from the default value in
#'   the pheatmap package, which can be found in the pheatmap documentation.
#'
#' @param breaks Optional sequence of numbers used to map heatmap values to
#'   colors.  Must cover the range from 0 to 1 and be one element longer than
#'   color vector. The default breaks are adapted to show relative abundance of
#'   major groups, while maintaining some contrast at the presence/absence
#'   level.  If value is NA (default for pheatmap), then the breaks are
#'   calculated automatically.
#'
#' @return A heatmap plot of the proportions of assigned lineages in each sample.
#'
#' @usage TaxonomicHeatmap <- function(otu.table, threshold=0, sample_order=NA,
#'  color=c("#FFFFFF", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE",
#'  "#08589E"), breaks=c(0, 0.00001, 0.001, 0.01, 0.10, 0.20, 0.30, 1), ...)
TaxonomicHeatmap <- function(otu.table, threshold=0, sample_order=NA,
                             color=c("#FFFFFF", "#CCEBC5", "#A8DDB5", "#7BCCC4", "#4EB3D3", "#2B8CBE", "#08589E"),
                             breaks=c(0, 0.00001, 0.001, 0.01, 0.10, 0.20, 0.30, 1),
                             ...) {
  assignment.counts <- LineageCounts(otu.table)

  # reorder the columns if requested
  if (!is.na(sample_order)) {
    assignment.counts <- assignment.counts[, sample_order]
  }
  
  # remove rows falling below the threshold
  if (threshold > 0) {
    assignment.sums <- apply(assignment.counts, 1, sum)
    assignment.counts <- assignment.counts[assignment.sums >= threshold,] 
  }

  assignment.fracs <- apply(assignment.counts, 2, function(x) { x / sum(x) })
  pheatmap(as.matrix(assignment.fracs), breaks=breaks, color=color, ...)
}


