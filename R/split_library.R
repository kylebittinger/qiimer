#' read split_library_log file and parse it.
#'
#' returns list with parsed quality metrics and Sample_Count dataframe
#' with SmapleID and read count for it.
#'
#' @param split_library_log_file log file
#' @param qiime_mapping_file map file
#' @note SampleID in mapping file and log file should be the same
#' @note stops if files do not exists
#' @note stops if cannot parse the table, or if SampleIDs are not the same 
read_split_library_log <- function(split_library_log_file, qiime_mapping_file) {
    if ( ! file.exists(split_library_log_file)) {
        stop("Cannot find split_library_log file.")
    }
    if ( ! file.exists(qiime_mapping_file)) {
        stop("Cannot find qiime_mapping_file file.")
    }

    totals <- parse_overall_quality(split_library_log_file)
    sample_read_count <- parse_sample_read_count_table(split_library_log_file)
    if (nrow(sample_read_count) == 0) {
        stop("Cannot find SampleID-Count table.")
    }
    map <- read_qiime_mapping_file(qiime_mapping_file)
    if ( ! is_the_same(sample_read_count$SampleID, map$SampleID)) {
        stop("Samples in split_library_log and qiime mapping file are different.")
    }

    c(totals, list("Sample_Count"=sample_read_count))
}

#' find lines that have given labels
parse_overall_quality <- function(split_library_log_file) {
    quality <- readLines(split_library_log_file)
    lines_labels <- c(
        "Total number of input sequences:",
        "Count of N characters exceeds limit:",
        "Barcode errors exceed max:",
        "Total number seqs written"
    )

    hits_against_labels <- sapply(lines_labels, 
        function(pattern) ifelse(grepl(pattern, quality), sub(pattern, "", quality), NA)
    )

    quality_metrics <- as.list(apply(hits_against_labels, 2, function(v) v[! is.na(v) ]))

    if (length(quality_metrics) != length(Filter(length, quality_metrics))) {
        stop("Cannot find all quality metrics.")
    }

    names(quality_metrics) <- c(
        "Total number of sequences",
        "Sequence errors exceed max",
        "Barcode errors exceed max",
        "Total number of sequences(after quality filtering)"
    )
    lapply(quality_metrics, function(m) as.numeric(m))
}

#' find SampleId - Read count info
parse_sample_read_count_table <- function(split_library_log_file, header_lines=15, tail_lines=2) {
    sample_count <- read.delim(split_library_log_file, header=FALSE, sep="\t", skip=header_lines)
    names(sample_count) <- c("SampleID", "ReadCount")
    sample_count$SampleID <- as.character(sample_count$SampleID)
    head(sample_count, -tail_lines)
}

is_the_same <- function(vec1, vec2) {
    setequal(vec1, vec2)
}
