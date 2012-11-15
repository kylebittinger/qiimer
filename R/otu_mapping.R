parse_otu_mapping <- function(filepath, prefix="") {
  # solution from http://stackoverflow.com/questions/6602881
  # Read in the data to character vector
  x <- scan(filepath, what="", sep="\n")
  # Separate elements by one or more whitepace
  y <- strsplit(x, "[[:space:]]+")
  # Extract the first vector element and set it as the list element name
  names(y) <- sapply(y, `[[`, 1)
  #names(y) <- sapply(y, function(x) x[[1]]) # same as above
  names(y) <- paste(prefix, names(y), sep="")
  # Remove the first vector element from each list element
  y <- lapply(y, `[`, -1)
  #y <- lapply(y, function(x) x[-1]) # same as above
  y
}

make_otu_table <- function(otus, sample_ids=NULL) {
  sample_vec <- factor(sub("_\\d+$", "", unlist(otus), perl=T))
  otu_vec <- factor(rep(names(otus), sapply(otus, length)))
  if (!is.null(sample_ids)) {
    extra_levels <- setdiff(levels(sample_vec), sample_ids)
    if (length(extra_levels) > 0) {
      extra_idx <- which(sample_vec %in% extra_levels)
      sample_vec <- factor(sample_vec[-extra_idx], levels=sample_ids)
      otu_vec <- factor(otu_vec[-extra_idx])
    } else {
      sample_vec <- factor(sample_vec, levels=sample_ids)
    }
  }
  table(otu_vec, sample_vec)
}
