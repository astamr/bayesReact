#' Function to process species-specififc sequence data and generate the seqs and seqList objects
#' @description Function using a path to a fasta file or seqs data frame to generate the seqs and seqList objects needed by Regmex to evaluate motif counts and sequence-specific probabilities, and saves the objects as independent .rds files. To obtain sequences from bioMart, please see the tutorial provided as part of the miReact package (https://github.com/muhligs/miReact/blob/master/makingUtrsForAdditionalSpecies.pdf).
#'
#' @param seq_in character or data frame, path to a bioMart output-file (decompressed .tsv file) or a data frame containing gene_name/gene_ID, sequences and sequence length.
#' @param out_path character, path to save seqs and seqList objects as .rds files. Default is current working directory.
#' @param gene_id character specifying the gene ID, which should be used (and this HAS to match the IDs in the gene/transcript expression matrix!). The input can be "gsym" (gene symbol, default), "gid" (ensembl gene id), "tid" (ensembl transcript id). Only relevant to specify when using a bioMart file-path as input.
#' @param min_length integer, minimum sequence length to include in the seqs object. Default is 20.
#' @param max_length integer, maximum sequence length to include in the seqs object. Default is 10,000.
#' @param one_seq_per_gene_id logical, if TRUE (default), only the longest sequence per gene will be included in the seqs object if a FASTA file is provided. Setting this parameter to FALSE is relevant when considering multiple transcripts per gene or when performing manual filtering of sequences.
#'
#' @importFrom dplyr %>%
#' @return Saves seqs and seqList objects as .rds files in the specified out_path.
#' @export
#'
#' @examples
#' \dontrun{
#' build_seq_list(seq_in = "./mart_export.txt", out_path = "./path/", gene_id = "gsym")
#' }
#'
build_seq_list <- function(seq_in, out_path = "./", gene_id = "gsym", min_length = 20, max_length = 10000, one_seq_per_gene_id = T){
  ## If bioMart output is provided as a path, load and process it ##
  if (is.character(seq_in)) {
    # initial checks (for suggested packages, which are not necessarily pre-installed)
    #if (isTRUE(suppressMessages(suppressWarnings(require(Biostrings)))) | isTRUE(suppressMessages(suppressWarnings(require(stringr)))))
    if (requireNamespace("Biostrings", quietly = TRUE) == F | requireNamespace("stringr", quietly = TRUE) == F) {
      stop("Please ensure to have both the 'Biostrings' & 'stringr' packages installed when using the function build_seq_list() with a fasta file as input.")
    }

    # load in the fasta file
    seqs_biostr <- suppressWarnings(Biostrings::readDNAStringSet(seq_in))
    seqs <- as.data.frame(seqs_biostr)

    # make columns from row name
    seqs[,2:5] <- stringr::str_split_fixed(rownames(seqs), "\\|", 4)
    seqs <- seqs[,c(2:5,1)]
    colnames(seqs) <- c("gid","tid","gsym", "strand", "sequence")
    seqs <- seqs[,c(gene_id, "sequence")]
    colnames(seqs) <- c("gid", "sequence")

    # Process and filter sequences
    seqs$sequence <- toupper(seqs$sequence) # sequence to upper case.
    seqs$nchar <- nchar(seqs$sequence) # sequence lengths

    # Some sequences have 'SEQUENCE' dummy annotation, which should be removed
    if (sum(!grepl("S",seqs$sequence)) > 0) {
      seqs$sequence <- sub("^SEQUENCE","",seqs$sequence) # make the dummy annotation have length 0
      seqs$nchar <- nchar(seqs$sequence) # re-define sequence lengths (will now be removed when conducting sequence length filtering).
    }
    # Sequences with Ns will interfere with downstream processes and should be removed
    if (sum(grepl("N",seqs$sequence)) > 0) {
      seqs <- seqs[!grepl("N",seqs$sequence),]
    }
    # Only retain one sequence per gene
    if (one_seq_per_gene_id) {
      if (requireNamespace("dplyr", quietly = TRUE) == F) {
        stop("Please ensure to have the 'dplyr' package installed when using the function build_seq_list() with a fasta file as input and one_seq_per_gene_id = TRUE.")
      }
      seqs <- seqs %>% group_by(gid) %>%
        filter(nchar == max(nchar)) %>%
        slice(1) %>% ungroup()
    }
    # Remove sequences with extremely short and long lengths (default is min_length < 20 and  max_length > 10000)
    seqs <- seqs[seqs$nchar > min_length & seqs$nchar < max_length,]
    print(paste0("After filtering, ", dim(seqs)[[1]], " sequences retained when constructing seqs data frame."))
    # Order by gene ID and length
    seqs <- seqs[order(seqs$gid,-seqs$nchar),]

    # Construct and save seq_list
    seqlist <- Regmex::seq.list.con(seqlist = seqs$sequence, cores=1)

    # Return seqs and seqlist
    return(list(seqs = seqs, seqlist = seqlist))

  ## If seqs DF is already generated, generate seqList ##
  } else if (is.data.frame(seq_in)) {
    # Check if seq_in contains correct columns
    if(!identical(colnames(seq_in), c("gid", "sequence", "nchar"))) stop("Input sequence data-frame must contain columns 'gid', 'sequence', and 'nchar'.", call. = F)
    # Filter sequences
    seqs <- seq_in[seq_in$nchar > min_length & seq_in$nchar < max_length,]
    # Generate seqlist
    seqlist <- Regmex::seq.list.con(seqlist = seq_in$sequence, cores=1)
    # Return seqs and seqlist
    return(list(seqs = seqs, seqlist = seqlist))

  } else { # check that this is only run if the other parts of the code is not run!
    stop("Sequence input must be either a path to a bioMart output-file or a data frame containing gene_name/gene_ID, sequences, and sequence lengths.")
  }
}


