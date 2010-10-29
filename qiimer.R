ReadMappingFile <- function(filepath) {
  mapping.file <- file(filepath)
  # Automatically closed during the call to read.table.

  # A "#" character may appear at start of header
  header <- readLines(mapping.file, n=1)
  header <- sub("^#", "", header)
  cols <- unlist(strsplit(header, "\t"))
  t <- read.table(mapping.file, col.names=cols, sep="\t", quote="",
                  row.names=NULL,
                  na.strings=c("NA", "na", "Null", "null"))

  t$SampleID <- as.character(t$SampleID)
  t$Description <- as.character(t$Description)
  t
}

ReadOtuTable <- function(filepath) {
  otu.file <- file(filepath)

  # Often, header line is preceeded by a comment line.
  headers <- readLines(otu.file, n=2)
  if (grep("Full OTU Counts", headers[1])) {
    header <- headers[2]
  } else {
    header <- headers[1]
  }
  cols <- unlist(strsplit(header, "\t"))

  cols[1] <- "OtuID"
  cols[length(cols)] <- "ConsensusLineage"

  t <- read.table(otu.file, col.names=cols, sep="\t", quote="",
                  row.names=NULL)

  t$OtuID <- as.character(t$OtuID)
  t
}
