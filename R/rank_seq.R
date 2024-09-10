#' Rank sequences based on fold-change (FC) score
#' @description Function to rank sequences based on FC-score, which is calculated as the difference between the expression of a sequence/gene and the median expression across all samples/cells.
#'
#' @param exp Numerical gene expression (or other relevant measure) matrix (Rows = genes, columns = samples). Input should NOT be log-transformed. Examples include raw counts (default), CPM, or TPM values.
#' @param data_type character, type of data to normalize and/or scale. Input can be "count" (default), "CPM" (counts per million; TPM can also be used), or "norm_scale_exp" (data already processed using norm_scale_seq()).
#' @param path character, directory to save .rds file in. Default is NULL, returning the fold-change score matrix.
#'
#' @return Sample/cell specific rank for each sequence, based on the fold-change (FC) score. Sequences with highest FC scores have the lowest ranks (thus the rank is decreasing).
#' @export
#'
#' @examples
#' \dontrun{
#' rank_of_seqs <- rank_seq(exp) # returns matrix with rank (integer) of each sequence for each sample.
#' }
#'
rank_seq <- function(exp, data_type = "count", path = NULL) {
  # Normalize and scale expression matrix
  if (data_type != "norm_scale_exp") {
    exp <- bayesReact::norm_scale_seq(exp, data_type = data_type, save_rds = F)
  }

  # Calculate FC scores
  medianExp <- apply(exp,1,stats::median)
  FC_rank <- apply(exp,2,function(x) order(x-medianExp, decreasing = TRUE))
  # Add rownames used to check overlap with seqs data.
  rownames(FC_rank) <- rownames(exp)

  if (is.character(path)) {
    file_path <- paste0(path, "FC_rank_", Sys.Date(), ".rds")
    saveRDS(FC_rank, file = file_path)
    return(file_path)
  }
  return(FC_rank)

  # !!!!!! Potentially setup so we either save FC_score or only compute and return medianExp! (depend on issues with data size)
}
