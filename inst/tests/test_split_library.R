context("read split_library_log file")

test_that("cannot read nonexistent file", {
    expect_error(
        read_split_library_log("fileNeverEverFoundAnyWhere.txt", "valid1_qiime_map.txt"), 
        "Cannot find split_library_log file.")
})

test_that("valid table has correct names", {
    qual <- read_split_library_log(
        "../testdata/valid1_split_library_log.txt",
        "../testdata/valid1_qiime_map.txt")
    expect_true(is.list(qual))
    expect_equal(names(qual), c(
        "Total number of sequences", 
        "Sequence errors exceed max", 
        "Barcode errors exceed max",
        "Total number of sequences(after quality filtering)",
        "Sample_Count")
    )
    expect_true(is.data.frame(qual[["Sample_Count"]]))
})

test_that("valid table has valid values for totals/errors/totals after split", {
    qual <- read_split_library_log(
        "../testdata/valid1_split_library_log.txt",
        "../testdata/valid1_qiime_map.txt")
    expect_equal(qual[["Sequence errors exceed max"]], 7250278)

    qual <- read_split_library_log(
        "../testdata/valid2_split_library_log.txt",
        "../testdata/valid2_qiime_map.txt")
    expect_equal(qual[["Total number of sequences"]], 597097)
    expect_equal(qual[["Sequence errors exceed max"]], 0)
    expect_equal(qual[["Barcode errors exceed max"]], 149480)
    expect_equal(qual[["Total number of sequences(after quality filtering)"]], 446734)
})

test_that("labels are in incorrect order but still can be read", {
    qual <- read_split_library_log(
        "../testdata/incorrect_order_split_library_log.txt",  
        "../testdata/valid1_qiime_map.txt") 
    expect_equal(qual[["Total number of sequences"]], 11190687)
    expect_equal(qual[["Sequence errors exceed max"]], 7250278)
    expect_equal(qual[["Barcode errors exceed max"]], 1110109)
    expect_equal(qual[["Total number of sequences(after quality filtering)"]], 2826739)
})

test_that("fails if some labels are not found", {
    expect_error(read_split_library_log(
        "../testdata/no_tot_no_barcode_split_library_log.txt", 
        "../testdata/valid1_qiime_map.txt"), 
        "Cannot find all quality metrics.")
})

test_that("can read SampleID-Read Count table", {
    qual <- read_split_library_log(
        "../testdata/valid1_split_library_log.txt",
        "../testdata/valid1_qiime_map.txt")
    expect_equal(ncol(qual$Sample_Count), 2)
    expect_equal(nrow(qual$Sample_Count), 65)

    qual <- read_split_library_log(
        "../testdata/valid2_split_library_log.txt",
        "../testdata/valid2_qiime_map.txt")
    expect_equal(ncol(qual$Sample_Count), 2)
    expect_equal(nrow(qual$Sample_Count), 16)
})

test_that("can read values of SampleID-Read Count table", {
    qual <- read_split_library_log(
        "../testdata/valid1_split_library_log.txt",
        "../testdata/valid1_qiime_map.txt")
    sample_count <- qual$Sample_Count
    has_read_count <- function(sample, count){
        hit <- data.frame(SampleID=sample, ReadCount=count)
        nrow(merge(hit, sample_count)) > 0
    }
    expect_true(has_read_count("PCRblank", 6))
    expect_true(has_read_count("GRO.49",32423))
    expect_true(has_read_count("Saliva2", 18))
})

test_that("fails if SampleID_Read Count table is not found", {
    expect_error(
        read_split_library_log(
            "../testdata/no_table_split_library_log.txt",
            "../testdata/valid1_qiime_map.txt"), 
        "Cannot find SampleID-Count table.")
})

test_that("fails if SampleID_Read Count table parsed has different Samples than mapping file", {
    expect_error(read_split_library_log(
        "../testdata/shifted_split_library_log.txt", 
        "../testdata/valid1_qiime_map.txt"), 
        "Samples in split_library_log and qiime mapping file are different.")
})

