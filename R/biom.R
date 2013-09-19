#' Extract biom data in raw form.
#' @param b An object of class `biom`.
#' @return For sparse biom objects, returns a data frame.
#'   For dense biom objects, returns a matrix.
#' @export
biom_raw_data <- function (b) {
  # Retrieve the first word from the biom type
  # One of: OTU, Taxon, Gene, Function, Ortholog, Pathway, Metabolite
  row_label <- strsplit(b$type, " ")[[1]][1]
  rnames <- sapply(b$rows, `[[`, "id")
  cnames <- sapply(b$columns, `[[`, "id")
  if (b$matrix_type == "sparse") {
    df <- as.data.frame(do.call(rbind, b$data))
    colnames(df) <- c(row_label, "SampleID", "Reads")
    df[[row_label]] <- factor(df[[row_label]], labels=rnames)
    df$SampleID <- factor(df$SampleID, labels=cnames)
    df
  } else if (b$matrix_type == "dense") {
    dnames <- list()
    dnames[[row_label]] <- rnames
    dnames[["SampleID"]] <- cnames
    matrix(unlist(b$data), nrow=b$shape[1], byrow=T, dimnames=dnames)
  } else {
    message("unsupported matrix type")
    NULL
  }
}

#' Extract taxonomy info from a biom object.
#' @param b An object of class `biom`.
#' @return A list of character vectors, one per row.
#' @export
biom_taxonomy <- function (b) {
  bt <- sapply(b$rows, function (x) x$metadata$taxonomy)
  names(bt) <- sapply(b$rows, `[[`, "id")
  bt
}