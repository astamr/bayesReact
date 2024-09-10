#' miReact code dependencies
#' @description This constitutes some of the core miReact code obtained from https://github.com/muhligs/miReact (19/12/2023).
#' For more information on miReact, please see their publication: https://doi.org/10.1038/s41598-021-88480-5 (Nielsen et al., Sci Rep, 2021).
#' bayesReact depends on both Regmex and miReact for handling input data and parallelization using the Slurm queuing system when running on a computer clusters.
#'
#' @param pattern state space (transition matrix) for motif of interest.
#' @param seq sequence for which the probability of motif occurrence is evaluated.
#'
#' @return motif probability (pd.mrs2)
#' or X (test2).
#' @export
#'
#' @examples
#' # See 'https://github.com/muhligs/miReact'.
#'
miReact_code <- function() NULL

#' @describeIn miReact_code Obtain motif probability in a given sequence.
pd.mrs2 <- function(pattern, seq){
  tm <- Regmex:::transition.matrix(pattern$matrix, seq$freq.mono)
  finl.st <- pattern$endState
  tm[finl.st,] <- 0
  tm[finl.st,finl.st] <- 1
  return(1-sum((tm %^% seq$length)[pattern$startState,-finl.st])) # %^% is from the expm package
}

#' @describeIn miReact_code Test.
test2 <- function(pattern, seq){
  2+2
}
