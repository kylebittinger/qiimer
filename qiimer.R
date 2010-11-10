library(reshape)

ReadCommentedHeader <- function(fileobj) {
  # Private function to read a (possibly) commented header line from a tsv file.
  
  header <- readLines(fileobj, n=1)
  header <- sub("^#", "", header)
  unlist(strsplit(header, "\t"))
}

ReadSampleMapping <- function(filepath) {
  # Read a QIIME sample mapping file and return a data frame.
  
  sample.file <- file(filepath, 'rt')

  cols <- ReadCommentedHeader(sample.file)

  sample.table <- read.table(sample.file,
                             col.names=cols,
                             sep="\t",
                             quote="",
                             row.names=NULL,
                             na.strings=c("NA", "na", "Null", "null"))
  close(sample.file)

  sample.table$SampleID <- as.character(sample.table$SampleID)
  if ('Description' %in% names(sample.table)) {
    sample.table$Description <- as.character(sample.table$Description)
  }

  sample.table
} 

GroupSamplesBy <- function(sample.mapping, colname) {
  tapply(sample.mapping$SampleID, sample.mapping[[colname]], c)
}

AttachMetadata <- function(otu.table, sample.mapping) {
  merge(otu.table, sample.mapping, by=c('SampleID'))
}

ReadOtuTable <- function(filepath) {
  # Read a QIIME OTU table and return a data frame.
  
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

# OTU table functions

SampleCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$SampleID, sum)
}
OtuCounts <- function(otu.table) {
  tapply(otu.table$Counts, otu.table$OtuID, sum)
}
LineageCounts <- function(otu.table) {
  tapply(otu.table$Counts, list(otu.table$Consensus.Lineage, otu.table$SampleID), sum)
}


ReadAssignmentTable <- function(filepath) {
  # Read a parsed table of OTU assignments

  read.table(filepath, header=TRUE, sep="\t", quote="", row.names=NULL)
}

kHeatmapBreaks <- c(0, 0.00001, 0.005, 0.01, 0.10, 0.20, 0.40, 0.60, 1)
kHeatmapColors <- c(rgb(255, 255, 255, max=255), rainbow(7))

TaxonomyHeatmap <- function(otu.table) {
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
            margins=c(5,25),
            cexRow=1,
            cexCol=1,
          )
  legend(0.8, 0.5,
         head(kHeatmapBreaks, n=-1),
         fill=kHeatmapColors,
         col=kHeatmapColors,
         )
}

