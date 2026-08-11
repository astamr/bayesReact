#' miReact code dependencies
#' @description This constitutes some of the core miReact code, which bayesReact depends on, obtained from https://github.com/muhligs/miReact (19/12/2023).
#' For more information on miReact, please see their publication: https://doi.org/10.1038/s41598-021-88480-5 (Nielsen et al., Sci Rep, 2021).
#' bayesReact depends on both Regmex and miReact for handling input data and parallelization using the Slurm queuing system when running on a computer clusters.
#'
#' @param pattern state space (transition matrix) for motif of interest.
#' @param seq sequence for which the probability of motif occurrence is evaluated.
#' @param markov_order Order of motif background model: `0` uses nucleotide
#'   frequencies (default), and `1` uses dinucleotide transition probabilities.
#'
#' @importFrom expm %^%
#' @return motif probability (pd_mrs2).
#'
#' @examples
#' # See 'https://github.com/muhligs/miReact'.
#'
#' @keywords internal
pd_mrs2 <- function(pattern, seq, markov_order = 0){ # in miReact, the function is named 'pd.mrs2'.
  if (markov_order == 0){
    #tm <- Regmex:::transition.matrix(pattern$matrix, seq$freq.mono)
    transition.matrix <- utils::getFromNamespace("transition.matrix", "Regmex")
    tm <- transition.matrix(pattern$matrix, seq$freq.mono)
    finl.st <- pattern$endState
    tm[finl.st,] <- 0
    tm[finl.st,finl.st] <- 1
    return(1-sum((tm %^% seq$length)[pattern$startState,-finl.st]))
  }
  if (markov_order == 1){
    #prob.dist.di <- utils::getFromNamespace("prob.dist.di", "Regmex")
    #return(prob.dist.di(pattern, seq, nt.null = 2, overlap = FALSE)$prob.1.or.more)
    transition.matrix.di <- utils::getFromNamespace("transition.matrix.di", "Regmex")

    tm <- transition.matrix.di(pattern$matrix.di, seq$con.prob.di)
    states <- rownames(tm)

    state_part <- sub("[ACGT]$", "", states)
    finl.st <- which(state_part %in% as.character(pattern$endState))

    tm[finl.st, ] <- 0
    tm[cbind(finl.st, finl.st)] <- 1

    init.st <- stats::setNames(numeric(length(states)), states)
    for (b in names(seq$freq.mono)) {
      di.st <- paste0(pattern$matrix[pattern$startState, b], b)
      init.st[di.st] <- init.st[di.st] + as.numeric(seq$freq.mono[b])
    }
    p_no_hit <- sum((init.st %*% (tm %^% (seq$length - 1L)))[, -finl.st, drop = FALSE])
    return(1 - p_no_hit)
  }
  stop("markov_order should either 0 or 1.", call. = F)
}

