#!/bin/bash
#SBATCH --account account_name
#SBATCH --mem 32G
#SBATCH -t 24:00:00

R -e "bayesReact::bayesReact_parallel(lst_data = list(FC_rank = \"/path/to/FC_rank.rds\",
                                                motif_probs = \"/path/to/motif_probs.rds\",
                                                motif_counts =  \"/path/to/motif_counts.rds\"),
                                      out_path = \"/path/to/out_folder/\",
                                      out_name = \"out_activity\",
                                      account = \"account_name\")"
