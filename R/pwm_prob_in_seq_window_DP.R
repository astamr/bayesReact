#' @keywords internal
#' @useDynLib bayesReact, .registration = TRUE
NULL

#' Fast sequence-specific PWM occurrence probability for single-sites (windows) using C++.
#' @description
#' Calls ./src/pwm_prob_dp.cpp
#'
#' Computes P(m), the probability that a PWM-width sequence window scores strictly above the specified cutoff,
#' given the nucleotide (markov_order = 0) or dinucleotide (markov_order = 1) composition of a sequence.
#'
#' @param pwm 4 x motif_width numeric matrix of log-odds scores, rows A, C, G, T.
#' @param cutoff numeric absolute score cutoff. A window is an occurrence
#' only when its summed score is strictly greater than this.
#' @param seqs_freq_nt n_seqs_eval x 4 numeric matrix of per-sequence nt frequencies.
#' @param cond_prob 4 x 4 x n_seqs_eval numeric array of dinucleotide conditional
#' probabilities, or NULL for a 0th-order background model.
#' @param resolution integer specifying number of score bins used to discretize the PWM score range.
#' Scores are rounded down during discretization, entailing a conservative estimate of P(m).
#' The approx. error in P(m) depends on how much probability mass are concentrated near the cutoff, and it decreases with increasing resolution.
#'
#' @keywords internal
pwm_prob_in_seq_window_DP <- function(pwm, cutoff, seqs_freq_nt, cond_prob = NULL, resolution = 10000) {
  # convert inputs to ones expected by C++
  storage.mode(pwm) <- "double"
  storage.mode(seqs_freq_nt) <- "double"
  if (!is.null(cond_prob)) storage.mode(cond_prob) <- "double"

  # Call C++ DP function
  window_prob <- .Call("pwm_prob_dp",
                       pwm, as.double(cutoff)[1L],
                       seqs_freq_nt, cond_prob,
                       as.integer(resolution)[1L],
                       PACKAGE = "bayesReact")

  return(window_prob) # P(m)
}
