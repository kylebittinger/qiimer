context('Import Sample Mapping file')

sample.mapping <- ReadSampleMapping('../testdata/sample_map.txt')

mapping.names <- c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Diet", "Description")
expect_that(names(sample.mapping), equals(mapping.names))

mapping.samples <- c("A.1", "A.2", "B.1", "B.2", "C.1")
expect_that(sample.mapping$SampleID, equals(mapping.samples))



context('Import OTU table')

otu.table <- ReadOtuTable('../testdata/otu_table.txt')

otu.ids <- as.factor(rep(seq(0, 4), 5))
expect_that(otu.table$OtuID, equals(otu.ids))

sample.ids <- as.factor(rep(mapping.samples, each=5))
expect_that(otu.table$SampleID, equals(sample.ids))




