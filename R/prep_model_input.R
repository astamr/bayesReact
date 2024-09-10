#' Helper function to prepare input data used for model fitting and motif activity estimation
#' @description
#' This function generate the appropriate input for the STAN model used to estimate motif activity.
#'
#' @param in_seq_motif_data list containing FC_rank (seq X sample), motif_probs (1 X seq), motif_counts (1 X seq).
#' @param threshold_motif_prob numeric value specifying the minimum threshold used to truncate the probability of observing a motif at least once in a sequence. If set to NULL, no threshold is used.
#' @param threshold_motif_count integer specifying the maximum threshold used to truncate the number of times a motif is observed in a sequence. If set to NULL, no threshold is used.
#'
#' @return List containing all necessary variables used as input for the STAN model.
#' @export
#'
#' @examples
#' \dontrun{
#' inlist <- prep_model_input(in_seq_motif_data =
#' list(FC_rank = FC_rank, motif_probs = motif_probs, motif_counts = motif_counts))
#' }
#'
prep_model_input <- function(in_seq_motif_data, threshold_motif_prob = 1e-10, threshold_motif_count = 2){
  # Data dimensions
  nr_obs <- dim(in_seq_motif_data$FC_rank)[[2]]
  nr_seqs <- dim(in_seq_motif_data$FC_rank)[[1]]

  # Re-scale total sequence-interval to have length one [0, 1] across which events/motifs occur
  if (is.numeric(threshold_motif_prob)) {
    l_vector <- (-log1p(-ifelse(in_seq_motif_data$motif_probs < threshold_motif_prob, threshold_motif_prob, in_seq_motif_data$motif_probs)))
  } else{
    l_vector <- (-log1p(-in_seq_motif_data$motif_probs)) # poisson lambda (depends on sequence specific probability of observing motif at least once)
  }

  l_vector <- l_vector/sum(l_vector) # scaled from 0-1

  l <- matrix(NA, nrow = nr_seqs, ncol = nr_obs)
  l[, 1:nr_obs] <- l_vector

  # Order seqs based on FC for each sample/cell
  l <- lapply(1:nr_obs, function(c) l[in_seq_motif_data$FC_rank[,c],c])
  l <- do.call(cbind, l)
  l <- rbind(rep(0, times = nr_obs), l) # add initial row due to indexing

  # cumsum for exact motif occurrence across the interval
  L <- apply(l, 2, function(i) cumsum(i))

  # generate indexing for motif occurrence
  s <- lapply(1:nr_obs, function(c) which(in_seq_motif_data$motif_counts[in_seq_motif_data$FC_rank[,c]] > 0))
  s <- do.call(cbind, s)

  n <- lapply(1:nr_obs, function(c) in_seq_motif_data$motif_counts[in_seq_motif_data$FC_rank[,c]][s[,c]])
  n <- do.call(cbind, n) # use to include intervals/seqs the number of times motif occurs

  s <- s+1 #  Indexing matching L

  # Create initial list
  d_list <- list(
    N = nr_seqs +1,
    K = dim(s)[[1]],
    C = nr_obs,
    l = l[,1:nr_obs],
    L = L[,1:nr_obs],
    s = s[,1:nr_obs],
    n = n[,1:nr_obs]
  )

  ########## generate model input data ##########

  r <- lapply(1:d_list$C, function(c) (d_list$L[d_list$s[,c],c] + d_list$L[d_list$s[,c]-1, c])/2)
  r <- do.call(cbind, r)

  ## pre-compute partial log-likelihood for each sample/cell for speed-up ##

  # precompute sum log length intervals
  l <- lapply(1:d_list$C, function(c) d_list$l[d_list$s[,c],c])
  l <- do.call(cbind, l)
  log_l <- d_list$n*log1p(l-1)
  sum_log_l <- colSums(log_l)

  # pre-compute partial log beta density (from likelihood)
  sum_log_r <- d_list$n*log1p(r-1)
  sum_log_r <- colSums(sum_log_r)
  sum_log_1_minus_r <- d_list$n*log1p(-r)
  sum_log_1_minus_r <- colSums(sum_log_1_minus_r)

  inlist <- list(
    K = sum(d_list$n[,1]),
    C = d_list$C,
    sum_log_l = sum_log_l,
    sum_log_r = sum_log_r,
    sum_log_1_minus_r = sum_log_1_minus_r
  )
  return(inlist)
}
