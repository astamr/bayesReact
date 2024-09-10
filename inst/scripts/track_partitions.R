#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

nr_jobs <- as.numeric(args[1]) #  total number of jobs
dir <- as.character(args[2]) # path to tmp work directory
out_path <- getwd() # path to output directory
setwd(paste0("./", dir)) # move to temp working directory
input_parameters <- readRDS("input_parameters.rds")
#logs_con <- file(input_parameters$logs_file)
logs_con <- input_parameters$logs_file
account <- input_parameters$account

## Track partitions ##
cat("\n___________Model fitting for data partitions____________\n\n", file = logs_con, append = T)
cat("All partitions have been submitted.\n", file = logs_con, append = T)
#cat("Jobs started: 0.\n", sep="", file = logs_con, append = T)
#cat("Jobs completed: 0.\n", sep="", file = logs_con, append = T)
cat("Jobs started: ", 0, ", Jobs completed: ", 0, sep="", file = logs_con, append = T) #"\n",
conn <- file(logs_con, open = "r+") # open log file in new connection
all_lines <- readLines(conn, warn = F)
close(conn)
last_line<- length(all_lines)
all_lines <- paste(all_lines[1:(last_line-1)], collapse = "\n")
all_lines <- paste0(all_lines, "\n")
#close(logs_con)
# monitor results until all are completed
job_compl <- length(list.files(pattern = "res*"))  # jobs completed
job_start <- length(list.files(pattern = "[.]out$")) # jobs started
job_status <- c(0,0) # job inspector
while(job_compl < nr_jobs){ # while not all jobs are completed
  Sys.sleep(20) # reporting delay
  if(!identical(c(job_start, job_compl), job_status)){ # if something new happened
    job_status <- c(job_start, job_compl) # update job status
    updateLog <- paste0("Jobs started: ", job_start, ", Jobs completed: ", job_compl) # Update job status
    cat(all_lines, file = input_parameters$logs_file) # rewrite log-file
    cat(updateLog, file = input_parameters$logs_file, append = T) # update last line
  }
  job_compl <- length(list.files(pattern = "res*"))  # update jobs completed
  job_start <- length(list.files(pattern = "[.]out$")) # update jobs started
}
# Write final progress out
updateLog <- paste0("Jobs started: ", job_start, ", Jobs completed: ", job_compl, ".\n")
cat(all_lines, file = input_parameters$logs_file) # rewrite log-file
cat(updateLog, file = input_parameters$logs_file, append = T) # update last line

## Start partition collection and clean-up job once all sub-jobs are finished ##

# Start job to collect all partitions
#logs_con <- file(input_parameters$logs_file)
cat("\nStarting job to collect results from all partitions.", "\n\n", sep="", file = logs_con, append = T)

# prepare SBATCH script
setwd(out_path)
sink("bayesCollect.script")
cat("#!/bin/sh", paste0("#SBATCH --time ", 12, ":00:00"), paste0("#SBATCH --mem ", 128, "gb"), paste0("#SBATCH --account ", account), "#SBATCH --job-name bayesCollect", sep="\n")
cat("\n[ $SLURM_PROCID -ne 0 ] && exit 0\n")
cat("cd ", out_path, "\n") #getwd()
cat(paste0("Rscript --vanilla ", base::system.file("", package = "bayesReact"), "//scripts/part_collect_and_clean_up.R "), out_path, " ", dir, "\n", sep="")
sink()
# run script
waste <- system("sbatch bayesCollect.script", ignore.stdout = T, ignore.stderr = T)
