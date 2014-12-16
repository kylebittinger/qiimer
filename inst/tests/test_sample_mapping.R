context('read_qiime_mapping_file')

s <- read_qiime_mapping_file(
  system.file("testdata", "sample_map.txt", package="qiimer"))

test_that("Column names are parsed correctly", {
  expect_equal(colnames(s), c(
    "SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Diet", 
    "Description"))
})

test_that("Sample IDs are parsed correctly", {
  expect_equal(class(s$SampleID), "character")
  expect_equal(s$SampleID, c("A.1", "A.2", "B.1", "B.2", "C.1"))
})

test_that("stop if file does not exist", {
    expect_error(read_qiime_mapping_file("NeverExistingFile.txt"), 
        "File NeverExistingFile.txt does not exist. Check the name or location of the file." );
})

test_that("stop if file has duplicate SampleIDs", {
    filepath <- system.file("testdata", "sample_map_duplicate_ids.txt", package="qiimer")
    expect_error(read_qiime_mapping_file(filepath), 
        paste("File",  filepath, "has duplicate SampleIDs. Check that SampleIDs are unique."));
})

test_that("Description column can be empty.", {
    filepath <- system.file("testdata", "sample_map_empty_description.txt", package="qiimer")
    s <- read_qiime_mapping_file(filepath);
    expect_true(all(is.na(s$Description)))
})

test_that("Non Description column can not be empty(for factors).", {
    filepath <- system.file("testdata", "sample_map_empty_non_description_column.txt", package="qiimer")
    expect_error(read_qiime_mapping_file(filepath),
        paste("Column(s)", "Diet", "have empty entry. Write explicitly NA or value needed."), 
        fixed=TRUE);
})

test_that("Non Description column can not be empty(for numbers).", {
    filepath <- system.file("testdata", "sample_map_empty_non_description_column_numeric.txt", package="qiimer")
    expect_error(read_qiime_mapping_file(filepath),
        paste("Column(s)", "Age", "have empty entry. Write explicitly NA or value needed."),
        fixed=TRUE)
})

