#' Read a QIIME OTU table and return a data frame.
ReadOtuTable <- function(filepath) {
  otu.file <- file(filepath, 'rt')

  # Often, header line is preceeded by a comment line.
  header <- readLines(otu.file, n=1)
  if (grepl("Full OTU Counts", header)) {
    header <- readLines(otu.file, n=1)
  }
  cols <- unlist(strsplit(header, "\t"))
  cols[1] <- "OtuID"
  
  otu.table <- read.table(otu.file,
                          col.names=cols,
                          sep="\t",
                          quote="",
                          as.is=TRUE,
                          row.names=NULL)
  close(otu.file)

  otu.table$OtuID <- as.factor(otu.table$OtuID)
  otu.table <- melt(otu.table,
                    id.vars=c('OtuID', 'Consensus.Lineage'),
                    variable_name="SampleID",
                    )
  otu.table <- rename(otu.table, c(value="Counts"))
  otu.table
}

#' Return a list of counts per sample.
SampleCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$SampleID, sum)
}

#' Return a list of counts per OTU.
OtuCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$OtuID, sum)
}

#' Return a list of counts per assigned lineage.
LineageCounts <- function(otu.table) {
  tapply(otu.table$Counts, list(otu.table$Consensus.Lineage, otu.table$SampleID), sum)
}

#' Attach sample metadata to OTU table data frame.
AttachMetadata <- function(otu.table, sample.mapping) {
  merge(otu.table, sample.mapping, by=c('SampleID'))
}

#' Default breakpoints for taxonomy heatmap.
kHeatmapBreaks <- c(0, 0.00001, 0.001, 0.01, 0.10, 0.20, 0.30, 1)

#' Default colors for taxonomy heatmap.
kHeatmapColors <- c(rgb(240, 249, 232, max=255),
                    rgb(204, 235, 197, max=255),
                    rgb(168, 221, 181, max=255),
                    rgb(123, 204, 196, max=255),
                    rgb(78, 179, 211, max=255),
                    rgb(43, 140, 190, max=255),
                    rgb(8, 88, 158, max=255))

#' Create a heatmap of taxonomic assignments.
#'
#' Additional arguments are passed directly to the pheatmap function.
TaxonomicHeatmap <- function(otu.table, breaks=kHeatmapBreaks, col=kHeatmapColors, ...) {
  assignment.counts <- LineageCounts(otu.table)
  assignment.fracs <- apply(assignment.counts, 2, function(x) { x / sum(x) })
  pheatmap(as.matrix(assignment.fracs), breaks=breaks, col=col, ...)
}

#' Create a histogram of proportional OTU counts.
#'
#' This can be useful for setting the breakpoints in an OTU Heatmap.
OtuHistogram <- function(otu.table) {
  assignment.counts <- LineageCounts(otu.table)
  assignment.fracs <- apply(assignment.counts, 2, function(x) { log10(x / sum(x)) })
  hist(assignment.fracs)
}

