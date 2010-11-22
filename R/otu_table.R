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

#' Return a list of counts per sample
SampleCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$SampleID, sum)
}

#' Return a list of counts per OTU
OtuCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$OtuID, sum)
}

#' Return a list of counts per assigned lineage
LineageCounts <- function(otu.table) {
  tapply(otu.table$Counts, list(otu.table$Consensus.Lineage, otu.table$SampleID), sum)
}

#' Attach sample metadata to OTU table data frame
AttachMetadata <- function(otu.table, sample.mapping) {
  merge(otu.table, sample.mapping, by=c('SampleID'))
}

#' Default breakpoints for taxonomy heatmap
kHeatmapBreaks <- c(0, 0.00001, 0.005, 0.01, 0.10, 0.20, 0.40, 0.60, 1)

#' Default colors for taxonomy heatmap
kHeatmapColors <- c(rgb(255, 255, 255, max=255), rainbow(7))

#' Create a heatmap of taxonomic assignments
TaxonomyHeatmap <- function(otu.table, cexRow=0.7) {
  assignment.counts <- LineageCounts(otu.table)
  assignment.fracs <- apply(assignment.counts, 2, function(x) { x / sum(x) })
  
  heatmap(as.matrix(assignment.fracs),

            # dendrogram control
            Rowv=TRUE,
            Colv=TRUE,
            
            scale="none",
            #labRow=assignment.labels,
            breaks=kHeatmapBreaks,
            col=kHeatmapColors,
            margins=c(5,18),
            cexRow=cexRow,
            cexCol=1,
          )

  legend(0.83, 0.5,
         head(kHeatmapBreaks, n=-1),
         fill=kHeatmapColors,
         col=kHeatmapColors,
         cex=0.6,
         )
}


