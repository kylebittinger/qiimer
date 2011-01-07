#' Read a QIIME sample mapping file.
#'
#' @param filepath Path to sample mapping file.  The file must conform to the
#'   QIIME standards, detailed at
#'   \url{http://qiime.org/documentation/file_formats.html}.
#'
#' @return A data frame.
ReadSampleMapping <- function(filepath) {  
  sample.file <- file(filepath, 'rt')

  header <- readLines(sample.file, n=1)
  header <- sub("^#", "", header)
  cols <- unlist(strsplit(header, "\t"))

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

  sample.table[order(sample.table$SampleID),]
} 

#' Group samples by common value in metadata column.
#'
#' @param sample.mapping Sample mapping data frame.
#'
#' @param colname Column name by which to group.
#'
#' @return A list of sample names per category.
GroupSamplesBy <- function(sample.mapping, colname) {
  tapply(sample.mapping$SampleID, sample.mapping[[colname]], c)
}

