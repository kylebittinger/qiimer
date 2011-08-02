#' Read a collated alpha diversity table from QIIME.
#'
#' @param filepath Path to alpha diversity table file.
#'
#' @return A data frame representing the table.
ReadAlphaDiversity <- function(filepath) {
  adiv.table <- read.table(filepath,
                           sep="\t",
                           quote="",
                           na.strings=c('NA', 'n/a'),
                           comment.char="#",
                           header=TRUE)
  adiv.table <- adiv.table[2:length(colnames(adiv.table))]
  colnames(adiv.table)[1] <- "num.seqs"
  adiv.table <- melt(adiv.table, id.vars=c("num.seqs", "iteration"), variable_name="SampleID")
  colnames(adiv.table)[length(colnames(adiv.table))] <- "diversity"
  adiv.table
}

#' Compute summary statistics for collated alpha diversity tables.
#'
#' @param adiv.table A collated alpha diversity data frame.
#'
#' @return A data frame of summary statistics.
AlphaDiversityStats <- function(adiv.table) {
  ddply(adiv.table,
        c("num.seqs", "SampleID"),
        summarise,
        mean=mean(diversity),
        sd=sd(diversity))
}
