#' Read a QIIME sample mapping file.
#'
#' @param filepath Path to sample mapping file.  The file must conform to the
#'   QIIME standards, detailed at
#'   \url{http://qiime.org/documentation/file_formats.html}.
#'
#' @return A data frame.
parse_mapping_file <- function(filepath) {  
  sample_file <- file(filepath, 'rt')

  header <- readLines(sample_file, n=1)
  header <- sub("^#", "", header)
  cols <- unlist(strsplit(header, "\t"))

  sample_data <- read.table(sample_file,
                             col.names=cols,
                             sep="\t",
                             quote="",
                             row.names=NULL,
                             na.strings=c("NA", "na", "Null", "null"))
  close(sample_file)

  sample_data$SampleID <- as.character(sample_data$SampleID)
  if ('Description' %in% names(sample_data)) {
    sample_data$Description <- as.character(sample_data$Description)
  }
  sample_data
} 


#' Group samples by common value in metadata column.
#'
#' @param sample.mapping Sample mapping data frame.
#'
#' @param colname Column name by which to group.
#'
#' @return A list of sample names per category.
group_samples_by_column <- function(sample.mapping, colname) {
  tapply(sample.mapping$SampleID, sample.mapping[[colname]], c)
}

