library('RUnit')
source('qiimer.R')

mapping.fp <- "tests/sample_map.txt"
otu.table.fp <- "tests/otu_table.txt"

test.suite <- defineTestSuite("qiimer",
                              dirs = file.path("tests"),
                              testFileRegexp = '^test.+\\.R')
test.result <- runTestSuite(test.suite)
printTextProtocol(test.result, showDetails=TRUE)
