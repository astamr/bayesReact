#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

out_path <- as.character(args[1]) # path to output directory
dir <- as.character(args[2]) # path to tmp work directory
system("rm trackPart*")
setwd(paste0("./", dir)) # move to temp working directory
input_parameters <- readRDS("input_parameters.rds")
output_type <- input_parameters$output_type # output type
model <- input_parameters$model # model
out_name <- input_parameters$out_name # output name
save_as_bigmat <- input_parameters$save_as_bigmat # save as bigmat
#logs_con <- file(input_parameters$logs_file)
logs_con <- input_parameters$logs_file

## Collecting partitions and generating output ##
cat("\n__________Collecting partitions and clean-up____________\n\n", file = logs_con, append = T)
part_res <- list.files(pattern = "res_")
cat("Loading ", length(part_res), " result files. \n", file = logs_con, append = T)
part_res <- part_res[order(as.numeric(sub("-.*", "", sub("res_", "", part_res))), decreasing = FALSE)] # define order for loading results

if (output_type == "activity_summary"){
  for(i in seq_along(part_res)){ # load results into motif_activities
    base::load(part_res[i])
    if(i==1){
      motif_a_mean <- motif_act$motif_a_mean
      motif_a_sd <- motif_act$motif_a_sd
      motif_activity <- motif_act$motif_activity
      motif_post_prob <- motif_act$motif_post_prob
      motif_a_CI_lower <- motif_act$motif_a_CI_lower
      motif_a_CI_upper <- motif_act$motif_a_CI_upper
      motif_model_nEff <- motif_act$motif_model_nEff
      motif_model_Rhat <- motif_act$motif_model_Rhat
      if(model == "BF"){
        motif_model_lBF <- motif_act$motif_model_lBF
      }
    }
    if (i>1){
      motif_a_mean <- rbind(motif_a_mean, motif_act$motif_a_mean)
      motif_a_sd <- rbind(motif_a_sd, motif_act$motif_a_sd)
      motif_activity  <- rbind(motif_activity, motif_act$motif_activity)
      motif_post_prob <- rbind(motif_post_prob, motif_act$motif_post_prob)
      motif_a_CI_lower <- rbind(motif_a_CI_lower, motif_act$motif_a_CI_lower)
      motif_a_CI_upper <- rbind(motif_a_CI_upper, motif_act$motif_a_CI_upper)
      motif_model_nEff <- rbind(motif_model_nEff, motif_act$motif_model_nEff)
      motif_model_Rhat <- rbind(motif_model_Rhat, motif_act$motif_model_Rhat)
      if(model == "BF"){
        motif_model_lBF <- rbind(motif_model_lBF, motif_act$motif_model_lBF)
      }
    }
  }
}

if (output_type == "activity"){
  for(i in seq_along(part_res)){ # load results into motif_activities
    base::load(part_res[i])
    if(i==1){
      motif_activities <- motif_act$motif_activity
    }
    if (i>1){
      motif_activities  <- rbind(motif_activities, motif_act$motif_activity)
    }
  }
}

## Cleaning up temporary files ##
cat("Cleaning up temporary files. \n", file = logs_con, append = T)
system("rm res_*")
#system("rm Rscript*")
system("rm bayesPart*")
system("rm sbatch.script")
system("rm *.rds")
setwd(out_path) # jump to old working directory (output path)
system(paste0("rmdir ", dir)) # remove temp working directory

## Save output ##
cat("Saving output. \n", file = logs_con, append = T)
if(save_as_bigmat == F){
  if(model == "BF"){
    motif_activities <- list(motif_activity = motif_activity, motif_post_prob = motif_post_prob, motif_a_mean = motif_a_mean, motif_a_sd = motif_a_sd,
                             motif_a_CI_lower = motif_a_CI_lower, motif_a_CI_upper = motif_a_CI_upper,
                             motif_model_nEff = motif_model_nEff, motif_model_Rhat = motif_model_Rhat, motif_model_lBF = motif_model_lBF)}
  if (output_type == "activity_summary" & model != "BF"){
    motif_activities <- list(motif_activity = motif_activity, motif_post_prob = motif_post_prob, motif_a_mean = motif_a_mean, motif_a_sd = motif_a_sd,
                             motif_a_CI_lower = motif_a_CI_lower, motif_a_CI_upper = motif_a_CI_upper,
                             motif_model_nEff = motif_model_nEff, motif_model_Rhat = motif_model_Rhat)
  }
  saveRDS(motif_activities, paste0(out_path, "/", out_name, ".rds")) # gzip compression by default
  cat(paste0("Output has succesfully been saved to ", out_path, "/", out_name, ".rds \n"), file = logs_con, append = T)
  #close(logs_con)
  #return(paste0(out_path, out_name, ".rds"))
  cat("\nDone. \n", file = logs_con, append = T)
} else{
  system(paste0("mkdir ", out_name))
  out <- paste0(out_path, "/", out_name, "/")
  # PROBABLY NEED TO START EARLIER IN BUILDING AND APPENDING BIG MATRICES IN ORDER TO OPTIMISE MEMORY USAGE!!!

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Should be file backed setup instead!!!! loads less into memory!!!!!
  # OR check if reading txt loads more into memory than .bin and .desc!!!??
  # Probably also better to write to file(s) iteratively
  if(output_type == "activity"){
    bigmemory::write.big.matrix(motif_activities, filename = paste0(out, "motif_activities.txt"), col.names = TRUE, row.names = TRUE)
  } else{
    bigmemory::write.big.matrix(motif_activity, filename = paste0(out, "motif_activity.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_post_prob, filename = paste0(out, "motif_post_prob.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_a_mean, filename = paste0(out, "motif_a_mean.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_a_sd, filename = paste0(out, "motif_a_sd.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_a_CI_lower, filename = paste0(out, "motif_a_CI_lower.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_a_CI_upper, filename = paste0(out, "motif_a_CI_upper.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_model_nEff, filename = paste0(out, "motif_model_nEff.txt"), col.names = TRUE, row.names = TRUE)
    bigmemory::write.big.matrix(motif_model_Rhat, filename = paste0(out, "motif_model_Rhat.txt"), col.names = TRUE, row.names = TRUE)
    if(model == "BF"){
      bigmemory::write.big.matrix(motif_model_lBF, filename = paste0(out, "motif_model_lBF.txt"), col.names = TRUE, row.names = TRUE)
    }
  }
  cat(paste0("Output has succesfully been saved to ", out_path, "/", out_name, " folder. \n"), file = logs_con, append = T)
  #close(logs_con)
  #return(paste0(out_path, "/", out_name, "/"))
  cat("\nDone. \n", file = logs_con, append = T)
}

# WORKING WITH BIGMEMORY MATRICES:
#motif_probs <- bigmemory::as.big.matrix(motif_probs, type = "double", backingfile = "motif_probs.bin", backingpath = paste0(getwd(), "/", dir), descriptorfile = "motif_probs.desc")
#motif_counts <- bigmemory::as.big.matrix(motif_counts, type = "double", backingfile = "motif_counts.bin", backingpath = paste0(getwd(), "/", dir), descriptorfile = "motif_counts.desc")
#cat("Writing motif probabilities and counts to disk. \n")
#motif_probs <- bigmemory::as.big.matrix(motif_probs, type = "double")
#motif_counts <- bigmemory::as.big.matrix(motif_counts, type = "double")
#bigmemory::write.big.matrix(motif_probs, filename = paste0(getwd(), "/", dir, "/motif_probs.txt"), col.names = TRUE, row.names = TRUE)
#bigmemory::write.big.matrix(motif_counts, filename = paste0(getwd(), "/", dir, "/motif_counts.txt"), col.names = TRUE, row.names = TRUE)

#bigmemory::write.big.matrix(motif_probs, filename = paste0(infile, "_motif_probs.txt"), col.names = TRUE, row.names = TRUE)
#bigmemory::write.big.matrix(motif_counts, filename = paste0(infile, "_motif_counts.txt"), col.names = TRUE, row.names = TRUE)
#system(paste0("cp ", "motif_probs.txt ", "./", infile, "_motif_probs.txt"))
#system(paste0("cp ", "motif_counts.txt ", "./", infile, "_motif_counts.txt"))

#motif_probs <- bigmemory::attach.big.matrix("motif_probs.desc", header = T, has.row.names = T, shared = T)
#motif_counts <- bigmemory::attach.big.matrix("motif_counts.desc", header = T, has.row.names = T, shared = T)
#motif_probs <- bigmemory::read.big.matrix(filename = paste0(partition, "_motif_probs.txt"), type = "double", header = TRUE, has.row.names = TRUE)
#motif_counts <- bigmemory::read.big.matrix(filename = paste0(partition, "_motif_counts.txt"), type = "double", header = TRUE, has.row.names = TRUE)
#motif_probs <- bigmemory::read.big.matrix(filename = "motif_probs.txt", type = "double", header = TRUE, has.row.names = TRUE)
#motif_counts <- bigmemory::read.big.matrix(filename = "motif_counts.txt", type = "double", header = TRUE, has.row.names = TRUE)
