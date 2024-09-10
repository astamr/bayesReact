#' Normalize and scale raw count matrix
#' @description Normalize total sample/cell expression to 1, multiply by median, and log2 transform.
#'
#' @param exp matrix, raw expression matrix (Rows = genes, columns = samples). CPM matrix can also be used, however, ensure that the data is NOT on the log-scale (?).
#' @param data_type character, type of data to normalize and/or scale. Input can be "count" (default) or "CPM" (counts per million; TPM can also be used).
#' @param save_rds logical, whether to save the normalized and scaled matrix as an .rds file. Default is TRUE.
#' @param path character, directory to save .rds file in. Default is current working directory.
#'
#' @return Normalized and scaled expression matrix (Rows = genes, columns = samples).
#' @export
#'
#' @examples
#' \dontrun{
#' norm_scale_exp <- norm_scale_seq(exp, save_rds = F)
#' # Returns norm and scaled exp matrix.
#'
#' norm_scale_exp_path <- norm_scale_seq(exp, save_rds = T)
#' # Saves norm matrix and returns path to .rds file.
#' }
#'
norm_scale_seq <- function(exp, data_type = "count", save_rds = T, path = "./") {
  # Initial check(s)
  if(NA %in% exp) stop("Input expression matrix contains NA values. Please remove rows containing NAs.", call. = F)
  if(!(data_type %in% c("count", "CPM", "TPM"))) stop("Input data_type is not valid. Please use either 'count', 'CPM', or 'TPM'.", call. = F)

  # Remove genes not expressed
  rm <- which(rowSums(exp) == 0)
  if (length(rm) != 0){
    exp <- exp[-rm,]
  }

  # Normalize and log2 transform
  if (data_type == "count"){
    # normalize use total gene count in sample as approx. library size, and multiply by median expression of each sample
    exp <- sweep(exp, 2, colSums(exp), FUN = '/')*stats::median(colSums(exp)) # normalize columns to 1 and multiply by median
    exp <- log2(exp+1)
  }
  if (data_type %in% c("CPM", "TPM")){ # TEST!!!!!!!!!!!!!!!!!
    # Scale expression by multiplying with the median of each sample
    exp <- apply(exp, 2, function(x) x*stats::median(x))
    exp <- log2(exp+1)
  }

  # Save or return normalized expression matrix
  if(save_rds) {
    file_path <- paste0(path, "norm_scale_exp_", Sys.Date(), ".rds")
    saveRDS(exp, file = file_path)
    return(file_path)
  }
  return(exp)
}
