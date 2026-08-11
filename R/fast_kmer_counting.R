#' @keywords internal
#' @useDynLib bayesReact, .registration = TRUE
NULL

#' Fast exact non-overlapping counting of all k-mers using C++.
#' @description
#' Calls ./src/fast_kmers_counting.cpp
#'
#' @keywords internal
count_all_kmers_nonoverlap <- function(k, seqlist) {
  k <- as.integer(k)

  # check valid input
  if (length(k) != 1 || is.na(k) || k < 1 || k > 15) {
    stop("'k' must be an integer from 1 to 15.", call. = FALSE)}

  # extract nucleotide sequences
  seqs <- vapply(seqlist, `[[`, character(1L), "sequence")

  # count all k-mer occurrences within each sequence
  counts <- .Call("fast_kmers_counting",
    seqs, k, PACKAGE = "bayesReact")

  colnames(counts) <- Regmex::all.mers(k)
  return(counts)
}
