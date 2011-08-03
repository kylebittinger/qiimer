context('Import Sample Mapping file')

sample_table <- parse_mapping_file('../testdata/sample_map.txt')
sample_names <- c("A.1", "A.2", "B.1", "B.2", "C.1")

test_that("Sample column names are parsed correctly", 
  expect_that(colnames(sample_table), equals(
    c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Diet",
      "Description"))))
test_that("Sample IDs are parsed correctly", 
  expect_that(sample_table$SampleID, equals(
    sample_names)))


context('Import OTU table')

otu_table <- parse_otu_table('../testdata/otu_table.txt')

test_that("OTU IDs are parsed correctly",
  expect_that(otu_table$otu_ids, equals(
    c("0", "1", "2", "3", "5"))))
test_that("Sample IDs are parsed correctly",
  expect_that(otu_table$sample_ids, equals(
    sample_names)))




