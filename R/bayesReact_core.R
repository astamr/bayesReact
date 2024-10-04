#' bayesReact method for motif activity inference
#' @description
#' Function for predicting motif activity from ranked sequence data using a simple Bayesian model implemented in STAN.
#' bayesReact_core can be run locally and used for motif activity inference for a set of motifs across a small data set.
#' For larger data sets, please consider bayesReact_parallel().
#'
#' @param lst_data list containing file paths or data frames, e.g.: list(FC_rank = "./FC_rank_date.rds", motif_probs = "./seqXmot_probs.rds", motif_counts = "./seqXmot_counts.rds").
#' @param threshold_motif_prob numeric value specifying the minimum threshold used to truncate the probability of observing a motif at least once in a sequence. If set to NULL, no threshold is used.
#' @param threshold_motif_count integer specifying the maximum threshold used to truncate the number of times a motif is observed in a sequence. If set to NULL, no threshold is used.
#' @param model A character string specifying the model to be used for inference:
#' Either "bayesReact" (default) or "BF" (Bayes Factor; when comparing beta model against the uniform null model).
#' @param output_type type of output to be returned; either "activity" (only outputs activity score); "activity_summary" (default; outputs activity score, posterior mean and sd for the underlying activity parameter 'a', log P(|a| <= 0 | data), credible intervals, and simple model diagnostics);
#' "full_posterior" (outputs all MCMC iterations after warm-up period); or "full_model" (returns full stanfit model object). "full_posterior" and "full_model" can only be obtained for a single motif at a time.
#' @param CI credible interval (CI) to be returned for activity estimates. Default is 80% CI: c(0.10, 0.90).
#' @param MCMC_iterations number of iterations to be run by the MCMC sampler.
#' @param MCMC_chains number of independent MCMC chains to be used.
#' @param MCMC_warmup initial iterations to be discarded for each chain as warm-up/burn-in.
#' @param MCMC_cores number of cores, default is equal to the number of chains (which is the maximum number of cores that can be utilized by STAN's MCMC sampler).
#' Alternatively consider parallel::detectCores().
#' @param MCMC_keep_warmup whether to keep the warm-up iterations or not.
#' @param parallel this parameter should never be changed manually and is used internally by bayesReact_parallel().
#'
#' @return Motif activity estimates in a format specified by output_type.
#' @export
#'
#' @examples
#' \dontrun{
#' motif_activities <- bayesReact_core(lst_data =
#' list(FC_rank = FC_rank, motif_probs = motif_probs, motif_counts = motif_counts))
#' }
#'
bayesReact_core <- function(lst_data, threshold_motif_prob = 1e-10, threshold_motif_count = 2, #!!!!! Consider being able to specify sampler (e.g. MCMC/HMC, vb, pigeons) !!!!!
                            model = "bayesReact", output_type = "activity_summary", CI = c(0.10, 0.90), #CI = c(0.005, 0.995),
                            MCMC_iterations = 3000, MCMC_chains = 3, MCMC_warmup = 500, MCMC_cores = MCMC_chains,
                            MCMC_keep_warmup = F, parallel = F) {
  # start logo
  if(!parallel){
    cat("\n")
    cat("Initiating ...\n")
    system(paste0("cat ", base::system.file("", package = "bayesReact"), "//art/ascii_logo.txt"))
    cat("                                core \n")
    cat("\n")
    cat("________________________________________________________\n________________________________________________________\n")
    cat("\n")
  } else{
    nr_motif_part <- lst_data$nr_motif_part # in parallel, track number of motif partitions
    motif_names <- lst_data$motif_names # iterate over motif names
    partition <- 1} # track current motif partition

  # check correct input type
  if (!is.list(lst_data)) stop("lst_data must be a list containing file paths (or data frames) in the format provided by bayesReact::process_raw_input():
                                 list(FC_rank = \"./FC_rank_date.rds\", motif_probs = \"./seqXmot_probs.rds\", motif_counts = \"./seqXmot_counts.rds\")", call. = F)
  if (!(model %in% c("bayesReact", "BF", "bayesReact_2param", "bayesReact_shrinkage"))) stop("model must be either \"bayesReact\" or \"BF\"", call. = F)

  ## Read in relevant data ##
  if (is.character(lst_data$FC_rank)) {FC_rank <- readRDS(lst_data$FC_rank)} else {FC_rank <- lst_data$FC_rank}
  if (is.character(lst_data$motif_probs)) {motif_probs <- readRDS(lst_data$motif_probs)} else {motif_probs <- lst_data$motif_probs}
  if (is.character(lst_data$motif_counts)) {motif_counts <- readRDS(lst_data$motif_counts)} else {motif_counts <- lst_data$motif_counts}

  if (!is.vector(motif_probs) && output_type %in% c("full_posterior", "full_model") ||
      isFALSE(tryCatch(ncol(motif_probs) == 1)) && output_type %in% c("full_posterior", "full_model")) stop("output_type must be \"activity\" or \"activity_summary\" when providing multiple motifs", call. = F)

  ## Function to generate input data and fit model for motif m ##
  act_motif_m <- function(in_seq_motif_data, motif){
    invisible(gc(reset = TRUE, full = TRUE)) # !!!!!!!!!
    #Sys.sleep(1) # to see if it solves gc() issue

    # start message
    cat(paste0("Fitting model for \'", motif, "\'. "))
    # generate input data to evaluate activity of motif m across all samples
    inlist <- bayesReact::prep_model_input(in_seq_motif_data, threshold_motif_prob = threshold_motif_prob, threshold_motif_count = threshold_motif_count)

    # if running in BF mode
    if (model == "BF"){
      # fit model for each sample in parallel to obtain the sample/cell-specific logml and BF
      BF_out <- parallel::mclapply(1:inlist$C, function(c) bayesReact::fit_motif_model(input = list(C = 1, K = inlist$K, sum_log_l = inlist$sum_log_l[c], sum_log_r = inlist$sum_log_r[c], sum_log_1_minus_r = inlist$sum_log_1_minus_r[c]),
                                                                                       model = model_stan, model_type = "BF", output_type = output_type, CI = CI, iterations = MCMC_iterations, chains = 1, warmup = MCMC_warmup, cores = MCMC_cores,
                                                                                       keep_warmup = MCMC_keep_warmup), mc.cores = MCMC_cores-1)
      # combine results for all samples
      cat("Succesfull model fit. \n")

      # if running bayesReact_parallel(), check if new motif partition needs to be loaded
      try(if(parallel & partition < nr_motif_part & motif == colnames(motif_probs)[ncol(motif_probs)]){
        partition <<- partition + 1
        cat(paste0("Loading motif partition ", partition, ".\n"))
        # replace current motif probs and counts partition and load in next one
        motif_probs <<- readRDS(paste0(partition, "_motif_probs.rds")) # Updates global variables
        motif_counts <<- readRDS(paste0(partition, "_motif_counts.rds"))}, silent = T)

      if(output_type == "activity_summary") return(do.call(rbind, BF_out))
      if(output_type == "activity") return(unlist(BF_out))
    }
    # fit model for motif m
    m_out <- bayesReact::fit_motif_model(input = inlist, model = model_stan, output_type = output_type, CI = CI,
                                iterations = MCMC_iterations, chains = MCMC_chains, warmup = MCMC_warmup, cores = MCMC_cores,
                                keep_warmup = MCMC_keep_warmup)

    # if running bayesReact_parallel(), check if next motif partition needs to be loaded
    try(if(parallel & partition < nr_motif_part & motif == colnames(motif_probs)[ncol(motif_probs)]){
      partition <<- partition + 1
      cat(paste0("Loading motif partition ", partition, ".\n"))
      # replace current motif probs and counts partition and load in next one
      motif_probs <<- readRDS(paste0(partition, "_motif_probs.rds")) # Updates global variables
      motif_counts <<- readRDS(paste0(partition, "_motif_counts.rds"))}, silent = T)
    return(m_out)
  }

  if (is.vector(motif_probs) & !parallel | isTRUE(tryCatch(ncol(motif_probs) == 1)) & !parallel){ # When only one motif is provided; Enables returning full posterior or model object.
    if (model == "BF" & requireNamespace("bridgesampling", quietly = TRUE) == F) stop("The bridgesampling package is required for model = \"BF\". Please install bridgesampling and try again.", call. = F)
    if(!is.vector(motif_probs)){
      motif_probs <- motif_probs[,1]
      motif_counts <- motif_counts[,1]
    }
    if (length(intersect(names(motif_probs), rownames(FC_rank))) != length(rownames(FC_rank))) stop(paste0(length(rownames(FC_rank)) - length(intersect(names(motif_probs), rownames(FC_rank))), " genes/transcripts in the fold-change data do not have a matching entry in the provided motif probs/counts data. Please check if the same gene IDs/names are used?") , call. = F)
    # match cols of motif_probs and motif_counts to rows of FC_rank
    gene_set <- rownames(FC_rank)
    #FC_rank <- FC_rank[gene_set,]
    motif_probs <- motif_probs[gene_set]
    motif_counts <- motif_counts[gene_set]
    # check if rownames of FC_rank and names for motif_probs are identical
    if(!identical(rownames(FC_rank), names(motif_probs))) stop("rownames(FC_rank) and names(motif_probs) are not identical. Please check if the same gene IDs/names are used?" , call. = F)

    # Define STAN model
    model_stan <- bayesReact::construct_motif_model(model = model)
    # Fit model for motif
    act_motif <- act_motif_m(in_seq_motif_data = list(FC_rank = FC_rank, motif_probs = motif_probs, motif_counts = motif_counts), "motif")

    if (output_type == "activity_summary"){rownames(act_motif) <- colnames(FC_rank)}
    if (output_type == "activity"){names(act_motif) <- colnames(FC_rank)}
    if (output_type == "full_posterior"){colnames(act_motif$a) <- paste("a", colnames(FC_rank), sep = "_")}
    # return output
    return(act_motif)
  }

  if (!parallel){
    ## Match motif_probs and motif_counts to FC_rank ##
    if (model == "BF" & requireNamespace("bridgesampling", quietly = TRUE) == F) stop("The bridgesampling package is required for model = \"BF\". Please install bridgesampling and try again.", call. = F)
    # check if rows contain motifs instead of genes
    if (is.null(rownames(motif_probs))) {
      cat("\nWarning: No rownames are provided. Row numbers are used instead. \n")
      rownames(motif_probs) <- 1:nrow(motif_probs)
      rownames(motif_counts) <- 1:nrow(motif_counts)
    }
    check_correct_rows <- gsub("[AGCT]", "", rownames(motif_probs))
    if(sum(!grepl("[A-Za-z]", check_correct_rows)) == nrow(motif_probs) & !is.vector(motif_probs)){
      cat("\nWarning: motif_probs and motif_counts appear to be transposed. Transposing back... \n")
      motif_probs <- t(motif_probs)
      motif_counts <- t(motif_counts)
    }
    # check that motif_probs and motif_counts have column/motif names
    if (is.null(colnames(motif_probs))) {
      cat("\nWarning: No motif names (colnames) are provided. Column numbers are used instead. \n")
      colnames(motif_probs) <- 1:ncol(motif_probs)
      colnames(motif_counts) <- 1:ncol(motif_counts)
    }
    # Check for duplicated motif names
    if(T %in% duplicated(colnames(motif_probs))){
      cat("\nWarning: Duplicated motif names detected. Using row numbers as column names for motif_probs & motif_counts. \n\n", file = logs_con, append = T)
      colnames(motif_probs) <- 1:ncol(motif_probs)
      colnames(motif_counts) <- 1:ncol(motif_counts)
    }
    # check if motif_probs and motif_counts have same number of rows as FC_rank
    if (length(intersect(rownames(motif_probs), rownames(FC_rank))) != length(rownames(FC_rank))) stop(paste0(length(rownames(FC_rank)) - length(intersect(rownames(motif_probs), rownames(FC_rank))), " genes/transcripts in the fold-change data do not have a matching entry in the provided motif probs/counts data. Please check if the same gene IDs/names are used?") , call. = F)
    # match rows of motif_probs and motif_counts to rows of FC_rank
    gene_set <- rownames(FC_rank)
    #FC_rank <- FC_rank[gene_set,]
    motif_probs <- motif_probs[gene_set,,drop=F]
    motif_counts <- motif_counts[gene_set,,drop=F]
    # check if rownames of FC_rank and motif_probs are identical
    if(!identical(rownames(FC_rank), rownames(motif_probs))) stop("rownames(FC_rank) and rownames(motif_probs) are not identical. Please check if the same gene IDs/names are used?" , call. = F)

    motif_names <- colnames(motif_probs)
  }

  # Threshold for motif count
  if(is.integer(threshold_motif_count)){
    motif_counts[motif_counts > threshold_motif_count] <- threshold_motif_count
  }

  ## Estimate activity for each motif m ##
  # Define STAN model
  model_stan <- bayesReact::construct_motif_model(model = model)

  # Fit model for each motif m
  act_motif <- lapply(motif_names, function(m) act_motif_m(in_seq_motif_data = list(FC_rank = FC_rank, motif_probs = motif_probs[,m], motif_counts = motif_counts[,m]), motif = m))
  # Generate output matrices
  if (output_type == "activity_summary"){
    motif_a_mean <- do.call(cbind, lapply(act_motif, function(x) x$mean))
    motif_a_sd <- do.call(cbind, lapply(act_motif, function(x) x$sd))
    motif_activity <- do.call(cbind, lapply(act_motif, function(x) x$activity))
    motif_post_prob <- do.call(cbind, lapply(act_motif, function(x) x$post_prob))
    motif_a_CI_lower <- do.call(cbind, lapply(act_motif, function(x) x[,paste0(CI[1]*100, "%")]))
    motif_a_CI_upper <- do.call(cbind, lapply(act_motif, function(x) x[,paste0(CI[2]*100, "%")]))
    motif_model_nEff <- do.call(cbind, lapply(act_motif, function(x) x$n_eff))
    motif_model_Rhat <- do.call(cbind, lapply(act_motif, function(x) x$Rhat))
    if (model == "BF"){
      motif_model_lBF <- do.call(cbind, lapply(act_motif, function(x) x$lBF))
    }
    # Add row and col names
    rownames(motif_a_mean) <- colnames(FC_rank)
    colnames(motif_a_mean) <- motif_names
    rownames(motif_a_sd) <- colnames(FC_rank)
    colnames(motif_a_sd) <- motif_names
    rownames(motif_activity) <- colnames(FC_rank)
    colnames(motif_activity) <- motif_names
    rownames(motif_post_prob) <- colnames(FC_rank)
    colnames(motif_post_prob) <- motif_names
    rownames(motif_a_CI_lower) <- colnames(FC_rank)
    colnames(motif_a_CI_lower) <- motif_names
    rownames(motif_a_CI_upper) <- colnames(FC_rank)
    colnames(motif_a_CI_upper) <- motif_names
    rownames(motif_model_nEff) <- colnames(FC_rank)
    colnames(motif_model_nEff) <- motif_names
    rownames(motif_model_Rhat) <- colnames(FC_rank)
    colnames(motif_model_Rhat) <- motif_names
    if (model == "BF"){
      rownames(motif_model_lBF) <- colnames(FC_rank)
      colnames(motif_model_lBF) <- motif_names
    }
    # Generate output list with matrices
    if(model == "BF"){
      out_matrices <- list(motif_activity = motif_activity, motif_post_prob = motif_post_prob, motif_a_mean = motif_a_mean, motif_a_sd = motif_a_sd,
                           motif_a_CI_lower = motif_a_CI_lower, motif_a_CI_upper = motif_a_CI_upper,
                           motif_model_nEff = motif_model_nEff, motif_model_Rhat = motif_model_Rhat, motif_model_lBF = motif_model_lBF)
    } else {
      out_matrices <- list(motif_activity = motif_activity, motif_post_prob = motif_post_prob, motif_a_mean = motif_a_mean, motif_a_sd = motif_a_sd,
                           motif_a_CI_lower = motif_a_CI_lower, motif_a_CI_upper = motif_a_CI_upper,
                           motif_model_nEff = motif_model_nEff, motif_model_Rhat = motif_model_Rhat)
    }
  }

  if (output_type == "activity"){
    motif_activity <- do.call(cbind, act_motif)
    # Add row and col names
    rownames(motif_activity) <- colnames(FC_rank)
    colnames(motif_activity) <- motif_names
    # Generate output list with matrices
    out_matrices <- list(motif_activity = motif_activity)
  }

  # Return matrices
  cat("Done.\n")
  return(out_matrices)
}

