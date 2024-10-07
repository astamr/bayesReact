#' bayesReact parallelization for motif activity inference on large data sets
#' @description
#' This function is a wrapper for bayesReact_core(), allowing for motif activity inference on large single cell atlases.
#' bayesReact_parallel() performs data partitioning and currently utilizes the Slurm job scheduler, and is thus designed to be run on a computer cluster.
#'
#'
#' @param lst_data list containing file paths to input data or data frames, e.g.: list(FC_rank = "./FC_rank_date.rds", motif_probs = "./seqXmot_probs.rds", motif_counts = "./seqXmot_counts.rds").
#' To generate input data, please see the process_raw_input() function.
#' @param out_path full path to output directory (e.g., "./out/path/"). If the output directory does not exist initially, it will be created.
#' @param out_name name of output file. Default is "motif_activity".
#' @param save_as_bigmat boolean specifying whether to save output as big.matrix object (default is FALSE). For data sets larger than 100,000 cells/samples, this parameter is recommended to be TRUE.
#' @param account Slurm project account or individual user account used to track resources for parallel jobs.
#' @param samples_per_partition number of samples to be processed in each parallel job. Default is 20. Smaller values will result in more jobs that individually will run faster, however, more jobs needs to be allocated resources by the Slurm queuing system.
#' @param threshold_motif_prob numeric value specifying the minimum threshold used to truncate the probability of observing a motif at least once in a sequence. If set to NULL, no threshold is used.
#' @param threshold_motif_count integer specifying the maximum threshold used to truncate the number of times a motif is observed in a sequence. If set to NULL, no threshold is used.
#' @param model A character string specifying the model to be used for inference:
#' Either "bayesReact" (default) or "BF" (Bayes Factor; to compare the performance of the beta model against the uniform null model).
#' @param output_type type of output to be returned; either "activity" (only outputs activity score), or "activity_summary" (default; outputs activity score, posterior mean and sd for the underlying activity parameter 'a', log P(|a| <= 0 | data), credible intervals, and simple model diagnostics).
#' @param CI credible interval (CI) to be returned for activity estimates. Default is 80% CI: c(0.10, 0.90).
#' @param MCMC_iterations number of iterations to be run by the MCMC sampler.
#' @param MCMC_chains number of independent MCMC chains to be used.
#' @param MCMC_warmup initial iterations to be discarded for each chain as warm-up/burn-in.
#' @param MCMC_cores number of cores, default is equal to the number of chains (which is the maximum number of cores that can be utilized by STAN's MCMC sampler).
#' @param MCMC_keep_warmup whether to keep the warm-up iterations or not.
#' @param posterior_approx algorithm for approximating the posterior distribution: "MCMC" (default) or "Laplace" (faster but less accurate approximation).
#' Only works with model = "bayesReact" and output_type = "activity".
#'
#' @return The function returns the path to the output file, which is saved in the 'out_path' and named 'out_name'. If output_type = "activity", the output contains a matrix of dimension motifs X cells/samples with motif activities.
#' If output_type = "activity_summary", the output contains a list with the motif X cell/sample matrices: "motif_activity" (the motif activity estimates), "motif_post_prob" (the posterior probability of log P(|a| <= 0 | data) + log(2)),
#' "motif_a_mean" (posterior mean of 'a'), "motif_a_sd" (posterior standard deviation of 'a'),
#' "motif_a_CI_lower" (lower bound of credible interval (CI) for the underlying activity parameter 'a'), motif_a_CI_upper" (upper bound of CI),
#' "motif_model_nEff" (effective sample sizes; autocorrelation diagnostic), and "motif_model_Rhat" (R-hat convergence diagnostic).
#' If model = "BF", then the output also contains the "motif_model_lBF" matrix (log Bayes Factor (BF) values).
#' @export
#'
#' @examples
#' \dontrun{
#' results_out_path <- bayesReact_parallel(lst_data = list(FC_rank = "./FC_rank_date.rds",
#' motif_probs = "./seqXmot_probs.rds", motif_counts = "./seqXmot_counts.rds"),
#' out_path = "./", out_name = "hs_motif_activity", account = "my_account_id")
#' # Here, the output is saved in the current working directory as "hs_motif_activity.rds"
#' }
#'
bayesReact_parallel <- function(lst_data, out_path, out_name = "motif_activity", save_as_bigmat = F, account, samples_per_partition = 20,
                                threshold_motif_prob = 1e-10, threshold_motif_count = 2,
                                model = "bayesReact", output_type = "activity_summary", CI = c(0.10, 0.90), #CI = c(0.005, 0.995),
                                MCMC_iterations = 3000, MCMC_chains = 3, MCMC_warmup = 500,
                                MCMC_cores = MCMC_chains, MCMC_keep_warmup = F, posterior_approx = "MCMC"){

  logs_file <- paste0(out_path, "logs_", out_name, ".txt")
  logs_con <- logs_file #file(logs_file, open = "a")
  cat("\n", file = logs_con)
  cat("Initiating ...\n", file = logs_con, append = T)
  system(paste0("cat ", base::system.file("", package = "bayesReact"), "//art/ascii_logo.txt", " >> ", logs_file))
  cat("                            parallel \n", file = logs_con, append = T)
  cat("\n", file = logs_con, append = T)
  cat("________________________________________________________\n________________________________________________________\n", file = logs_con, append = T)
  cat("\n", file = logs_con, append = T)

  ## Checks ##
  if (!(model %in% c("bayesReact", "BF", "bayesReact_2param"))) {
    cat("Error: \'model\' must be either \"bayesReact\" or \"BF\"\n", file = logs_con, append = T)
    stop("model must be either \"bayesReact\" or \"BF\"", call. = F)
  }
  if (!(posterior_approx %in% c("MCMC", "Laplace"))) {
    cat("Error: \'posterior_approx\' must be either \"MCMC\" or \"Laplace\"\n", file = logs_con, append = T)
    stop("posterior_approx must be either \"MCMC\" or \"Laplace\"", call. = F)
  }
  if (posterior_approx == "Laplace" & model != "bayesReact" & output_type != "activity") {
    cat("Error: Laplace approximation only works with the \"bayesReact\" model specification and \"activity\" output\n", file = logs_con, append = T)
    stop("Laplace approximation only works with the \"bayesReact\" model specification and \"activity\" output", call. = F)
  }
  if (model == "BF" & requireNamespace("bridgesampling", quietly = TRUE) == F) {
    cat("Error: The bridgesampling package is required for model = \"BF\". Please install bridgesampling and try again.\n", file = logs_con, append = T)
    stop("The bridgesampling package is required for model = \"BF\". Please install bridgesampling and try again.", call. = F)
  }
  input_error_message <- "lst_data must be a list containing file paths in the format provided by bayesReact::process_raw_input():
                          list(FC_rank = \"./FC_rank_date.rds\", motif_probs = \"./seqXmot_probs.rds\", motif_counts = \"./seqXmot_counts.rds\")"
  if (!is.list(lst_data)) {
    cat(paste0("Error: ", input_error_message, "\n"), file = logs_con, append = T)
    stop(input_error_message, call. = F)
  }
  if (F %in% (c("FC_rank", "motif_probs", "motif_counts") %in% names(lst_data))) {
    cat(paste0("Error: ", input_error_message, "\n"), file = logs_con, append = T)
    stop(input_error_message, call. = F)
  }
  cat("\n___________Loading and partitioning input data__________\n\n", file = logs_con, append = T)
  if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = T)
    cat(paste0("Created output directory: ", out_path, " \n"), file = logs_con, append = T)
  }
  if (substring(out_path, nchar(out_path)) != "/"){out_path <- paste0(out_path, "/")}

  ## Read in relevant data ##
  input_parameters <- as.list(environment()) # to be saved and used by call_bayesReact_core()
  if (is.character(lst_data$FC_rank)) {FC_rank <- readRDS(lst_data$FC_rank)} else {FC_rank <- lst_data$FC_rank}
  if (is.character(lst_data$motif_probs)) {motif_probs <- readRDS(lst_data$motif_probs)} else {motif_probs <- lst_data$motif_probs}
  if (is.character(lst_data$motif_counts)) {motif_counts <- readRDS(lst_data$motif_counts)} else {motif_counts <- lst_data$motif_counts}
  if (output_type %in% c("full_posterior", "full_model")) {
    cat("Error: output_type must be \"activity\" or \"activity_summary\" when running bayesReact_parallel()", file = logs_con, append = T)
    stop("output_type must be \"activity\" or \"activity_summary\" when running bayesReact_parallel()", call. = F)
  }
  if (is.null(colnames(motif_probs))) {
    cat("\nWarning: No motif names (colnames) are provided. Column numbers are used instead. \n\n", file = logs_con, append = T)
    colnames(motif_probs) <- 1:ncol(motif_probs)
    colnames(motif_counts) <- 1:ncol(motif_counts)
  }

  ## Check if rows contain motifs instead of genes ##
  if (is.null(rownames(motif_probs))) {
    cat("\nWarning: No rownames are provided. Row numbers are used instead. \n\n", file = logs_con, append = T)
    rownames(motif_probs) <- 1:nrow(motif_probs)
    rownames(motif_counts) <- 1:nrow(motif_counts)
  }
  check_correct_rows <- gsub("[AGCT]", "", rownames(motif_probs))
  if(sum(!grepl("[A-Za-z]", check_correct_rows)) == nrow(motif_probs) & !is.vector(motif_probs)){
    cat("\nWarning: motif_probs and motif_counts appear to be transposed. Transposing back... \n\n", file = logs_con, append = T)
    motif_probs <- t(motif_probs)
    motif_counts <- t(motif_counts)
  }
  # Check for duplicated motif names
  if(T %in% duplicated(colnames(motif_probs))){
    cat("\nWarning: Duplicated motif names detected. Using row numbers as column names for motif_probs & motif_counts. \n\n", file = logs_con, append = T)
    colnames(motif_probs) <- 1:ncol(motif_probs)
    colnames(motif_counts) <- 1:ncol(motif_counts)
  }

  nr_obs <- dim(FC_rank)[[2]]
  if(nr_obs > 100000 & save_as_bigmat == F){
    cat("\nWarning: The number of observations is greater than 100,000. Consider setting save_as_bigmat = T to save the output as a memory efficient big.matrix object. \n\n", file = logs_con, append = T)
  }
  if(save_as_bigmat == T){
    if(requireNamespace("bigmemory", quietly = TRUE) == F) {
      cat("Error: The bigmemory package is required for save_as_bigmat = T. Please install bigmemory and try again.", file = logs_con, append = T)
      stop("The bigmemory package is required for save_as_bigmat = T. Please install bigmemory and try again.", call. = F)
    }
  }

  ## Ensure motif_probs and motif_counts are matrices (in case of single motif provided) ##
  if (is.vector(motif_probs)){
    motif_probs <- matrix(motif_probs, ncol = 1, dimnames = list(names(motif_probs), c("1")))
    motif_counts <- matrix(motif_counts, ncol = 1, dimnames = list(names(motif_counts), c("1")))
  }
  nr_motifs <- dim(motif_probs)[[2]]

  ## Match motif_probs and motif_counts to FC_rank ##
  # check if motif_probs and motif_counts have same number of rows as FC_rank
  if (length(intersect(rownames(motif_probs), rownames(FC_rank))) != length(rownames(FC_rank))) {
    cat(paste0("Error: ", length(rownames(FC_rank)) - length(intersect(rownames(motif_probs), rownames(FC_rank))), " genes/transcripts in the fold-change data do not have a matching entry in the provided motif probs/counts data. Please check if the same gene IDs/names are used?"), file = logs_con, append = T)
    stop(paste0(length(rownames(FC_rank)) - length(intersect(rownames(motif_probs), rownames(FC_rank))), " genes/transcripts in the fold-change data do not have a matching entry in the provided motif probs/counts data. Please check if the same gene IDs/names are used?") , call. = F)}
  # match cols of motif_probs and motif_counts to rows of FC_rank
  gene_set <- rownames(FC_rank)
  FC_rank <- FC_rank[gene_set,]
  motif_probs <- motif_probs[gene_set,, drop = F]
  motif_counts <- motif_counts[gene_set,, drop = F]
  # check if rownames of FC_rank and motif_probs are identical
  if(!identical(rownames(FC_rank), rownames(motif_probs))) {
    cat("Error: rownames(FC_rank) and rownames(motif_probs) are not identical. Please check if the same gene IDs/names are used?", file = logs_con, append = T)
    stop("rownames(FC_rank) and rownames(motif_probs) are not identical. Please check if the same gene IDs/names are used?", call. = F)
  }

  setwd(out_path) # set working directory as the output directory and create a temporary folder for handling temporary partition files and results.
  dir <- paste0("tmp_bayesReact_", out_name) # create temporary working directory
  #relaunch <- dir.exists(dir) # check if temporary working directory already exists
  if (dir.exists(dir)) {
    cat("\nWarning: Temporary working directory already exists and partial results present will be re-used. If different settings are supplied, please remove temporary folder before running bayesReact again! \n\n", file = logs_con, append = T)
    }else{
      cat("Creating temporary working directory ./", dir, "\n", sep="", file = logs_con, append = T)
      system(paste0("mkdir ", dir))}
  setwd(paste0("./", dir)) # move to temp working directory

  ## Make motif partitions ##
  motifs_per_partition <- min(nr_motifs, 500)
  part_start <- seq(1, nr_motifs, by = motifs_per_partition)
  part_end <- unique(c(seq(motifs_per_partition, nr_motifs, by = motifs_per_partition), nr_motifs))
  cat("Saving ", length(part_end),  " temporary motif partitions.", "\n", sep="", file = logs_con, append = T)
  for (i in seq_along(part_end)){
    saveRDS(motif_probs[,part_start[i]:part_end[i], drop = F], paste0(i, "_motif_probs.rds"))
    saveRDS(motif_counts[,part_start[i]:part_end[i], drop = F], paste0(i, "_motif_counts.rds"))
  }
  input_parameters$nr_motif_part <- length(part_end) # report number motif partitions to be loaded by bayesReact_core() one at a time (to reduce memory usage)
  input_parameters$motif_names <- colnames(motif_probs) # iterated over by bayesReact_core()
  saveRDS(input_parameters, "input_parameters.rds")

  ## Make job partitions ##
  # make intervals to partition samples (FC_rank) by
  samples_per_partition <- min(samples_per_partition, nr_obs)
  part_start <- seq(1, nr_obs, by = samples_per_partition)
  part_end <- unique(c(seq(samples_per_partition, nr_obs, by = samples_per_partition), nr_obs))
  tmp_file_names <- paste0(part_start, "-", part_end)
  # Specify job resources
  t <- 4 + round(8*(nr_motifs)/(4^7), 0) + as.integer(MCMC_iterations/10000)
  mem <- 3 + as.integer((MCMC_iterations*MCMC_chains)/10000) + as.integer(samples_per_partition/500)
  if (model == "BF") {
    #t <- 6 + round(8*(nr_motifs)/(4^7), 0) + as.integer(MCMC_iterations/10000)
    mem <- 4 + as.integer((MCMC_iterations*MCMC_chains)/10000) + as.integer(samples_per_partition/500)
  }
  # save partition files
  cat("Saving temporary data partitions and submitting ", length(part_end), " jobs.", "\n", sep="", file = logs_con, append = T)
  for(i in seq_along(part_end)){ # loop through the intervals for the FC_rank partitions
    infile <- paste0("FC_", tmp_file_names[i]) # generate file name for partition i
    # create partition
    if(!file.exists(infile)){ # check if partition file already exists, from previous partial run
      part_FC_rank <- FC_rank[,part_start[i]:part_end[i], drop = F]
      assign(infile, part_FC_rank)
      save(list=infile, file=paste0(infile,".Rdata")) # save in tmp wd for worker to load
    }

    ## Start job ##
    # prepare SBATCH script
    sink("sbatch.script")
    #cat("#!/bin/sh", "#SBATCH --nodes 1", paste0("#SBATCH -c ", MCMC_cores), "#SBATCH --time 12:00:00", "#SBATCH --mem 8gb", paste0("#SBATCH --account ", account), "#SBATCH --job-name Rscript", sep="\n")
    #cat("#!/bin/sh", "#SBATCH --nodes 1", paste0("#SBATCH -c ", MCMC_cores), "#SBATCH --time 12:00:00", "#SBATCH --mem-per-cpu 3gb", paste0("#SBATCH --account ", account), "#SBATCH --job-name Rscript", sep="\n")
    cat("#!/bin/sh", "#SBATCH --nodes 1", paste0("#SBATCH -c ", MCMC_cores), paste0("#SBATCH --time ", t, ":00:00"), paste0("#SBATCH --mem-per-cpu ", mem, "gb"), paste0("#SBATCH --account ", account), "#SBATCH --job-name bayesPart", sep="\n")
    cat("\n[ $SLURM_PROCID -ne 0 ] && exit 0\n")
    cat("cd ", getwd(), "\n")
    cat(paste0("Rscript --vanilla ", base::system.file("", package = "bayesReact"), "//scripts/call_bayesReact_core.R "), infile, " ", dir, "\n", sep="")
    sink()
    # run script
    waste <- system("sbatch sbatch.script", ignore.stdout = T, ignore.stderr = T)
  }

  ## Monitoring job progress ##
  # Start small tracking job
  cat("Initiating tracking of partition jobs.", "\n\n", sep="", file = logs_con, append = T)
  #close(logs_con)
  setwd(out_path)
  nr_jobs <- length(part_start) # total number of jobs
  sink("trackPart.script")
  cat("#!/bin/sh", paste0("#SBATCH -c ", 1), paste0("#SBATCH --time ", 48, ":00:00"), paste0("#SBATCH --mem ", 500, "MB"), paste0("#SBATCH --account ", account), "#SBATCH --job-name trackPart", sep="\n")
  cat("\n[ $SLURM_PROCID -ne 0 ] && exit 0\n")
  cat("cd ", out_path, "\n") #getwd()
  cat(paste0("Rscript --vanilla ", base::system.file("", package = "bayesReact"), "//scripts/track_partitions.R "), nr_jobs, " ", dir, "\n", sep="")
  sink()
  waste <- system("sbatch trackPart.script", ignore.stdout = T, ignore.stderr = T) # run tracking job

  return(invisible(NULL))
}
