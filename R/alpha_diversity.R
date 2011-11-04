#' Read a collated alpha diversity table from QIIME.
#'
#' @param filepath Path to alpha diversity table file.
#' @return A data frame representing the table.  The table is returned in
#'   melted form so sample IDs appear in rows rather than in columns.  The
#'   column of diversity values is named "diversity".
#' @export
parse_rarefaction <- function(filepath) {
  rarefaction_table <- read.table(filepath,
                           sep="\t",
                           quote="",
                           na.strings=c('NA', 'n/a'),
                           comment.char="#",
                           header=TRUE)
  rarefaction_table <- rarefaction_table[2:length(colnames(rarefaction_table))]
  colnames(rarefaction_table)[1] <- "sequences_per_sample"
  rarefaction_table <- melt(rarefaction_table, id.vars=c("sequences_per_sample", "iteration"), variable_name="SampleID")
  colnames(rarefaction_table)[length(colnames(rarefaction_table))] <- "diversity"
  rarefaction_table
}

#' Compute summary statistics for collated alpha diversity tables.
#'
#' @param rarefaction_table A collated alpha diversity data frame.
#' @return A data frame of summary statistics.
#' @export
rarefaction_stats <- function(rarefaction_table) {
  ddply(rarefaction_table, c("sequences_per_sample", "SampleID"),
        summarise, mean=mean(diversity), sd=sd(diversity))
}
