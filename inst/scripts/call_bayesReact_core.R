#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

partition <- as.character(args[1])
file_path <- paste0(partition, ".Rdata") # partition input file
dir <- as.character(args[2]) # path to work directory
#setwd(dir)
input_parameters <- readRDS("input_parameters.rds")
motif_probs <- readRDS("1_motif_probs.rds")
motif_counts <- readRDS("1_motif_counts.rds")

# set rstan options
rstan::rstan_options(auto_write = TRUE)

## Load FC partition and check if output is already produced ##
cat("Load input", format(Sys.time()), "\n")

if(file.exists(sub("FC_", "res_", file_path))){
  system(paste0("rm ", file_path))
  stop("Output already produced. Input removed.")
}
base::load(file_path) # loads part_FC_rank
FC_rank <- get(sub("[.]Rdata", "", file_path)) # redefine part_FC_rank as FC_rank
nr_obs <- dim(FC_rank)[[2]]

cat("Work-space size ", sum(sort( sapply(ls(), function(x) {object.size(get(x))})))/(1024*1024*1024), " GB ", date())
cat("\n")
cat("Running ", nr_obs, "samples", format(Sys.time()),"\n")

## Run bayesReact_core ##
motif_act <- bayesReact::bayesReact_core(lst_data = list(FC_rank = FC_rank, motif_probs = motif_probs, motif_counts = motif_counts,
                                                         nr_motif_part = input_parameters$nr_motif_part, motif_names = input_parameters$motif_names),
                                         threshold_motif_prob = input_parameters$threshold_motif_prob, threshold_motif_count = input_parameters$threshold_motif_count,
                                         model = input_parameters$model, parallel = T, output_type = input_parameters$output_type, CI = input_parameters$CI,
                                         MCMC_iterations = input_parameters$MCMC_iterations, MCMC_chains = input_parameters$MCMC_chains, MCMC_warmup = input_parameters$MCMC_warmup,
                                         MCMC_cores = input_parameters$MCMC_cores, MCMC_keep_warmup = input_parameters$MCMC_keep_warmup, posterior_approx = input_parameters$posterior_approx)

## Save partition output ##
cat("Saving results", sub("FC_", "res_", file_path), format(Sys.time()), " \n")
save(motif_act, file = sub("FC_","res_", file_path))
system(paste0("rm ", file_path))
cat(format(Sys.time()), "\n")
cat("Finished partition. \n")
