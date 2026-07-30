#' @keywords internal
#' @useDynLib bayesReact, .registration = TRUE
NULL

#' Fast exact non-overlapping counting of all k-mers using C++.
#'
#' @keywords internal
count_all_kmers_nonoverlap <- function(k, seqlist) {
  k <- as.integer(k)

  if (length(k) != 1L || is.na(k) || k < 1L || k > 15L) {
    stop("'k' must be an integer from 1 to 15.", call. = FALSE)
  }

  seqs <- vapply(seqlist, `[[`, character(1L), "sequence")

  counts <- .Call(
    "fast_kmers_counting",
    seqs,
    k,
    PACKAGE = "bayesReact")

  colnames(counts) <- Regmex::all.mers(k)
  counts
}
