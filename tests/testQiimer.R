test.ReadMappingFile <- function() {
  df <- ReadMappingFile(mapping.fp)
  checkEquals(names(df), c("SampleID", "BarcodeSequence",
                           "LinkerPrimerSequence", "Diet", "Description"))
  checkEquals(df$SampleID, c("A.1", "A.2", "B.1", "B.2", "C.1"))
}

test.ReadOtuTable <- function() {
  df <- ReadOtuTable(otu.table.fp)
  checkEquals(names(df), c("OtuID", "A.1", "A.2", "B.1", "B.2", "C.1",
                           "ConsensusLineage"))
  checkEquals(df$OtuID, c("0", "1", "2", "3", "4"))
}
