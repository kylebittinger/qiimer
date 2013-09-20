

otu_ids <- c("0", "1", "2", "3", "5")
sample_ids <- c("A.1", "A.2", "B.1", "B.2", "C.1")

if (suppressMessages(require(biom))) {

  context('biom_raw_data, dense format')

  db <- read_biom(
    system.file("testdata", "otu_table_dense.biom", package="qiimer"))

  test_that("OTU IDs are parsed correctly", {
    expect_equal(rownames(biom_raw_data(db)), otu_ids)
  })

  test_that("Sequence IDs are parsed correctly", {
    expect_equal(colnames(biom_raw_data(db)), sample_ids)
  })

  test_that("OTU counts are correct", {
    expect_equal(
      biom_raw_data(db)["0",], 
      structure(c(5, 10, 2, 0, 3), names=sample_ids))
    expect_equal(
      biom_raw_data(db)[,"A.1"], 
      structure(c(5, 76, 2, 637, 0), names=otu_ids))
  })

  context('biom_raw_data, sparse format')

  sb <- read_biom(
    system.file("testdata", "otu_table_sparse.biom", package="qiimer"))

  test_that("OTU IDs are present", {
    expect_equal(levels(biom_raw_data(sb)$OTU), otu_ids)
  })

  test_that("Sequence IDs are present", {
    expect_equal(levels(biom_raw_data(sb)$SampleID), sample_ids)
  })

  test_that("OTU counts are in the same order as the biom file", {
    expect_equal(
      biom_raw_data(sb)$Reads, 
      c(5, 10, 2, 3, 76, 23, 28, 43, 56, 2, 1, 637, 61, 63, 77, 34, 1))
  })

}
context('biom_taxonomy')

test_that("Taxonomy is correct", {
  expect_equal(biom_taxonomy(sb)[[1]], "Bacteria")
})

test_that("Taxonomy list is labeled with OTU IDs", {
  expect_equal(biom_taxonomy(sb)[["0"]], "Bacteria")
})
