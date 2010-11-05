ReadCommentedHeader <- function(fileobj) {
  # Private function to read a (possibly) commented header line from a tsv file.
  
  header <- readLines(fileobj, n=1)
  header <- sub("^#", "", header)
  unlist(strsplit(header, "\t"))
}

ReadSampleMapping <- function(filepath) {
  # Read a QIIME sample mapping file and return a data frame.
  
  mapping.file <- file(filepath)
  # Automatically closed during the call to read.table.

  cols <- ReadCommentedHeader(mapping.file)
  t <- read.table(mapping.file, col.names=cols, sep="\t", quote="",
                  row.names=NULL,
                  na.strings=c("NA", "na", "Null", "null"))

  t$SampleID <- as.character(t$SampleID)
  t$Description <- as.character(t$Description)
  t
}

ReadOtuTable <- function(filepath) {
  # Read a QIIME OTU table and return a data frame.
  
  otu.file <- file(filepath)

  # Often, header line is preceeded by a comment line.
  headers <- readLines(otu.file, n=2)
  if (grep("Full OTU Counts", headers[1])) {
    header <- headers[2]
  } else {
    header <- headers[1]
  }
  cols <- unlist(strsplit(header, "\t"))

  cols[1] <- "OTU"
  cols <- sub("Consensus Lineage", "Lineage", cols)
  
  t <- read.table(otu.file, col.names=cols, sep="\t", quote="",
                  row.names=NULL)

  t$OtuID <- as.character(t$OtuID)
  t
}

ReadAssignmentTable <- function(filepath) {
  # Read a parsed table of OTU assignments

  read.table(filepath, header=TRUE, sep="\t", quote="", row.names=NULL)
}

TaxonomyHeatmap <- function(otu.counts, otu.labels) {
  # Requires a modified OTU table containing taxa classified by rank.
  # This can be created by parse_taxonomy.py in the qiime_pipeline
  # module.
  #

  breaks <- c(-1, 0, 0.002, 0.004, 0.006, 0.01, 0.03, 0.05, 0.07, 0.09, 0.2, 0.4,
              0.6, 0.7, 0.9, 2)

  # set up color table
  colors <- c(rgb(255, 255, 255, max=255), rainbow(14))
  colors[7] <- rgb(0, 204, 0, max=255)
  colors[8] <- rgb(0, 255, 237, max=255)
  colors[11] <- rgb(25, 0, 217, max=255)
  colors[12] <- rgb(157, 0, 255, max=255)
  colors[13] <- rgb(220, 91, 255, max=255)
  colors[14] <- rgb(235, 182, 255, max=255)
  

  sample.labels <- names(otu.counts)
  assignments <- aggregate(otu.counts, list(Label=otu.labels), sum)

  assignment.labels <- assignments$Label
  assignment.counts <- assignments[,sample.labels]
  ConvertToProportion <- function(x) { x / sum(x) }
  assignment.fracs <- sapply(assignment.counts, ConvertToProportion)


  heatmap.2(as.matrix(assignment.fracs),

            # dendrogram control
            dendrogram="both",
            Rowv=TRUE,
            Colv=TRUE,
            
            scale="none",
            labRow=assignment.labels,
            breaks=breaks,
            col=colors,
            margins=c(5,20),
            cexRow=1,
            cexCol=1,
            density.info="none",
            trace="none",
          )

}
