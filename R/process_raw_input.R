#' Function to process raw input data needed to run bayesReact
#' @description
#' The function processes raw expression, motif and sequence data, in order to produce sequence ranks and motif probabilities needed for numerical sequence representation.
#' The processed data is used by bayesReact to evaluate motif distribution across the expression (or other continuous measure) ranked sequences.
#' The function calls the following sub-functions: norm_scale_seq, rank_seq, build_seq_list, and motif_prob, which can also be run individually.
#'
#' @param exp numerical gene expression (or other relevant measure) matrix (rows = genes, columns = samples) or character string specifying file path to expression data (saved in .rds format). Input should NOT be log-transformed. Examples include raw counts (default), CPM, or TPM values.
#' @param seq_data character string, data frame or vector; path to a bioMart output-file (decompressed .fasta file); a data frame containing gene_name/gene_ID, sequences and sequence length;
#' or a character vector of length two, containing the paths to seqs and seqlist objects (e.g., seq_data = c("./seqs.rds", "./seqlist.rds")).
#' @param motifs  vector of motifs (e.g., c("ACGTAGT", "GTACAAG")) or an integer specifying the k-mers to be used (e.g., 7 for all 7-mers).
#' A motif is specified as a regular expression on the alphabet {A,C,G,T}.
#' @param out_path character string specifying directory to save .rds files in. Default is current working directory.
#' @param exp_type character string specifying type of data to normalize and/or scale. Input can be "count" (default) or "CPM" (counts per million; TPM can also be used).
#' @param save_processed_exp logical, whether to also save the normalized and scaled matrix. Default is FALSE.
#' @param process_all logical, whether to process all data (expression, sequence, and motif data) or only computing the sample-specific fold-change ranks.
#' @param seq_gene_id character specifying the gene ID that should be used to build the sequence list. The IDs has to match the IDs in the gene/transcript expression matrix.
#' The input can be "gsym" (gene symbol, default), "gid" (ensembl gene id), "tid" (ensembl transcript id). Only relevant to specify when using a bioMart file-path as input.
#' @param seq_min_length integer, minimum sequence length to include in the seqs object and list. Default is 20.
#' @param seq_max_length integer, maximum sequence length to include in the seqs object and list. Default is 10,000.
#' @param approx_motif_prob logical, whether to use a binomial approximation to compute motif probabilities (default is FALSE).
#' @param cores integer, number of cores to use for parallel processing when computing motif probabilities (default is all available cores).
#'
#' @return process_raw_input() calculates and saves the Fold-change ranks, processed sequence data and sequence-specific motif probabilities, and returns a list containing all file paths; list(FC_rank_path, seqs_path, seqlist_path, motif_prob_path).
#' @export
#'
#' @examples
#' \dontrun{
#' out_file_paths <- process_raw_input(exp = "./project/folder/exp_counts.rds", seq_data = seqs,
#' motifs = 7, out_path = "./project/folder/")
#' }
#'
process_raw_input <- function(exp, seq_data, motifs, out_path = "./",
                              exp_type = "count", save_processed_exp = F, process_all = T,
                              seq_gene_id = "gsym", seq_min_length = 20, seq_max_length = 10000,
                              approx_motif_prob = F, cores = parallel::detectCores()){
  # initial checks
  if(is.data.frame(exp) + is.matrix(exp) + is.character(exp) == 0) stop("expression data input should be either a data frame, matrix, or character string specifying an .rds input file" , call. = F)

  ## Process expression data ##
  if (is.character(exp)){
    exp <- readRDS(exp)
  }
  # process expression data
  if (save_processed_exp){
    exp_path <- bayesReact::norm_scale_seq(exp, data_type = exp_type, path = out_path)
    exp <- readRDS(exp_path)
  } else {
    exp <- bayesReact::norm_scale_seq(exp, data_type = exp_type, save_rds = F)
  }

  if (process_all){
    ## Process sequence data ##
    if (!is.vector(seq_data)) {
      seq_out <- bayesReact::build_seq_list(seq_data, out_path = out_path, gene_id = seq_gene_id, min_length = seq_min_length, max_length = seq_max_length)
      seqs <- seq_out$seqs
      seqlist <- seq_out$seqlist

      # check overlap in sequences between expression and sequence data
      if (length(intersect(seqs$gid, rownames(exp)))/length(rownames(exp)) < 0.5) stop("Less than 50% of the genes in the expression data have a matching sequences in the provided sequence data. Please check if the same gene IDs/names are used?" , call. = F)
      cat(paste0(length(intersect(seqs$gid, rownames(exp)))/length(rownames(exp))*100, "% of genes in the expression data have a matching sequences in the provided sequence data. ",
                   length(rownames(exp) - length(intersect(seqs$gid, rownames(exp))), " genes will be removed from the expression data and fold-change calculation. \n")))

      # save seqs and seq list objects
      saveRDS(seqs, file = paste0(out_path, "seqs.rds"))
      saveRDS(seqlist, file = paste0(out_path, "seqlist.rds"))
      seq_paths <- list(seqs_path = paste0(out_path, "seqs.rds"), seqlist_path = paste0(out_path, "seqlist.rds"))
    }
    cat("Successfully processed sequence data. \n")

    ## Process motif data ##
    if (is.vector(seq_data)){
      seqs <- readRDS(seq_data[1])
      seqlist <- readRDS(seq_data[2])
      seq_paths <- list(seqs_path = seq_data[1], seqlist_path = seq_data[2])

      # check overlap in sequences between expression and sequence data
      if (length(intersect(seqs$gid, rownames(exp)))/length(rownames(exp)) < 0.5) stop("Less than 50% of the genes in the expression data have a matching sequences in the provided sequence data. Please check if the same gene IDs/names are used?" , call. = F)
      print(paste0(length(intersect(seqs$gid, rownames(exp)))/length(rownames(exp))*100, "% of genes in the expression data have a matching sequences in the provided sequence data. ",
                   length(rownames(exp) - length(intersect(seqs$gid, rownames(exp))), " genes will be removed from the expression data and fold-change calculation.")))

      # check that paths are given in correct order
      if (!is.data.frame(seqs)) stop("seq_data[1] input path should be to a data frame, please check that paths to sequence and sequence list objects haven't been swapped.", call. = F)
    }
    motif_prob_path <- bayesReact::motif_prob(motifs, seqs, seqlist, paths = F,
                                              cores = cores, out_path = out_path,
                                              binom_approx = approx_motif_prob, include_counts = !approx_motif_prob) # UPDATE include_counts IF FASTER APPROACH IS ADDED TO motif_prob()!!!

    cat("Successfully processed sequence-specific motif probabilities. \n")
  }
  gene_set <- intersect(seqs$gid, rownames(exp))
  exp <- exp[gene_set,]

  # generate fold-change based sequence ranks
  FC_rank_path <- list(FC_rank_path = bayesReact::rank_seq(exp, data_type = "norm_scale_exp", path = out_path))
  # If only processing expression data, return file path for the FC-based sequence ranks
  if (process_all == F) {
    warning("Only fold-change ranks were computed. Please set process_all = TRUE to process sequence and motif data. \nPlease ensure that you only included genes/transcripts for fold-change calculation, which are also present in the sequence data.")
    return(FC_rank_path)
  }
  cat("Successfully computed fold-change ranks. \n")

  # Save file paths for later use
  file_paths_out <- c(FC_rank_path, seq_paths, motif_prob_path)
  #saveRDS(file_paths_out, file = paste0(out_path, "data_paths.rds"))

  ## Return file paths ##
  return(file_paths_out)
}


