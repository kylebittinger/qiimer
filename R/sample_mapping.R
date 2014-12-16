#' Read a QIIME sample mapping file.
#'
#' @param filepath Path to sample mapping file.  The file must conform to the
#'   QIIME standards, detailed at
#'   \url{http://qiime.org/documentation/file_formats.html}.
#' @return A data frame of sample information.  Because the SampleID column is
#'   so often used to extract data from distance matrices and OTU tables, it 
#'   is returned as a character vector.
#' @export
read_qiime_mapping_file <- function(filepath) {  
  if ( ! file.exists(filepath)) {
    stop(paste("File", filepath, "does not exist. Check the name or location of the file."));
  }
  sample_file <- file(filepath, 'rt')

  header <- readLines(sample_file, n=1)
  header <- sub("^#", "", header)
  cols <- unlist(strsplit(header, "\t"))

  sample_data <- read.table(
    sample_file,
    col.names=cols,
    sep="\t",
    quote="",
    row.names=NULL,
    na.strings=c("NA", "na", "Null", "null"))
  close(sample_file)

  # The SampleID column is often used to extract data from distance matrices 
  # and OTU tables.  Storing as a character vector facilitates this practice.
  sample_data$SampleID <- as.character(sample_data$SampleID)
  if ( ! is_unique_ids(sample_data$SampleID)) {
    stop(paste("File", filepath, "has duplicate SampleIDs. Check that SampleIDs are unique."))
  }
  
  # Ideally, the Description column should contain a unique free text
  # description of each sample.  In this case, there is no advantage to using
  # a factor data type.
  if ('Description' %in% names(sample_data)) {
    sample_data$Description <- as.character(sample_data$Description)
  }

  invalid_columns <- column_with_empty_value(filepath, cols)
  if ( ! is.null(invalid_columns)) {
    stop(paste("Column(s)", invalid_columns, "have empty entry. Write explicitly NA or value needed."))
  }

  sample_data
} 

is_unique_ids <- function(data) {
    ids <- unlist(data)
    ids <- ids[! is.na(ids)]
    
    len <- length(ids)
    len_uniq <- length(unique(ids))

    (len == len_uniq)
}

# no implicit empty elements, but have to read file again
column_with_empty_value <- function(qiime_map_file, column_names, dont_check=c("Description")) {
    qiime_map <- read.table(
      qiime_map_file,
      col.names=column_names,
      sep="\t",
      quote="",
      row.names=NULL,
      colClasses="character",
      na.strings=c("NA", "na", "Null", "null")
    )
    columns_to_check <- setdiff(column_names, dont_check)
    is_valid_column <- sapply(columns_to_check, function(column) all(qiime_map[[column]] != ""))
    invalid_columns <- columns_to_check[ ! is_valid_column]
    if (length(invalid_columns) > 0) {
        return(invalid_columns)
    } 
    NULL
}
