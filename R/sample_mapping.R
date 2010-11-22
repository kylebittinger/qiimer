#' Read a possibly commented header line.
#'
ReadCommentedHeader <- function(fileobj) {
  header <- readLines(fileobj, n=1)
  header <- sub("^#", "", header)
  unlist(strsplit(header, "\t"))
}

#' Read a QIIME sample mapping file and return a data frame.
#'
#' @param filepath Path to sample mapping file.
ReadSampleMapping <- function(filepath) {  
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

#' Group samples by common value in metadata column.
#'
#' @param sample.mapping Sample mapping data frame.
#' @param colname Column name by which to group.
GroupSamplesBy <- function(sample.mapping, colname) {
  tapply(sample.mapping$SampleID, sample.mapping[[colname]], c)
}

