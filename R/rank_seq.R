#' Rank sequences based on fold-change (FC) score
#' @description Function to rank sequences based on FC-score, which is calculated as the difference between the expression of a sequence/gene and the median expression across all samples/cells.
#'
#' @param exp numerical gene expression (or other relevant measure) matrix (rows = genes, columns = samples). Input should NOT be log-transformed. Examples include raw counts (default), CPM, or RPKM values.
#' @param data_type character, type of data to normalize and/or scale. Input can be "count" (default), "CPM" (counts per million; TPM/RPKM can also be used), or "norm_scale_exp" (data already processed using norm_scale_seq()).
#' @param path character, directory to save .rds file in. Default is NULL, returning the fold-change score matrix.
#' @param control_exp numerical vector or matrix containing expression values for a user-defined control setting(s), matching the row ordering of the 'exp' input (vector) as well as the column ordering (matrix). If NULL (default), the median gene expressions are used instead.
#'
#' @return sample/cell specific rank for each sequence, based on the fold-change (FC) score. Sequences with highest FC scores have the lowest ranks (thus the rank is decreasing).
#' @export
#'
#' @examples
#' \dontrun{
#' rank_of_seqs <- rank_seq(exp) # returns matrix with rank (integer) of each sequence for each sample.
#' }
#'
rank_seq <- function(exp, data_type = "count", path = NULL, control_exp = NULL) {
  # normalize and scale expression matrix
  if (data_type != "norm_scale_exp") {
    exp <- bayesReact::norm_scale_seq(exp, data_type = data_type, save_rds = F)
  }

  # define the control setting
  if(is.null(control_exp)){
    medianExp <- apply(exp,1,stats::median)
  }else{
    # check that control_exp has the correct length
    if((is.vector(control_exp) && nrow(exp) != length(control_exp)) || (is.matrix(control_exp) && !identical(dim(exp), dim(control_exp)))) stop("provided expression profile(s) of the control setting(s) does not have the same dimensions as the 'exp' input" , call. = F)
    medianExp <- control_exp
  }

  # calculate FC scores
  if (is.vector(medianExp)) {
    FC_rank <- apply(exp,2,function(x) order(x-medianExp, decreasing = TRUE))
    # add rownames used to check overlap with seqs data
    rownames(FC_rank) <- rownames(exp)
  } else{
    # column-wise substraction of control condition
    FC_rank <- sapply(seq_len(ncol(exp)), function(c) {
      order(exp[, c] - medianExp[, c], decreasing = TRUE)
    })
    # add col & rownames
    colnames(FC_rank) <- colnames(exp)
    rownames(FC_rank) <- rownames(exp)
  }

  # return output
  if (is.character(path)) {
    if(substring(path, nchar(path)) != "/"){path <- paste0(path, "/")}
    file_path <- paste0(path, "FC_rank_", Sys.Date(), ".rds")
    saveRDS(FC_rank, file = file_path)
    return(file_path)
  }
  return(FC_rank)
}
