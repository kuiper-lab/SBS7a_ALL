#Load packages
source("~/git/pmc_kuiper_projects/HypermutatedALL/CODE/local/createdFunctions_HM-ALL.R")
library(tidyverse)
library(GenomicRanges)
library(MutationalPatterns)

# load rData
# load the merged mutation matricies of the cohort
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/SBS7aClusters_MutationMatrix.rdata")
# load the denovo extracted refitted signatures of the cohort
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/DeNovo_Refit_SBS7aALLClusters_ReuseSignatures.rdata")
#Load the list of SNVs
SNV_list <- read.table("~/Projects/hypermutated_ALL/ANALYSES/pathogenicVariants/20251217_clustered_SNV_list_pathogenic.tsv")

# for each extracted signature rescale the contribution to [0,1]
scaled_contributions <- reshape2::melt(fit_res_strict$contribution)
colnames(scaled_contributions) = c("Signature", "Sample", "Contribution")
samples <- unique(scaled_contributions$Sample)
scaled_contributions$Relative_contribution <- NA
for (i in 1:length(samples)){
  sample_specific_slice <- scaled_contributions[scaled_contributions$Sample == samples[i], ]
  sample_specific_slice_total_contribution <- sum(sample_specific_slice$Contribution)
  relative_contribution_vector <- sample_specific_slice$Contribution/sample_specific_slice_total_contribution
  scaled_contributions$Relative_contribution[scaled_contributions$Sample == samples[i]] <- relative_contribution_vector
}

# reshape the matrix to be better readable: we have the original and the rescaled
# values next to eachtother
scaled_contributions_df <- reshape(scaled_contributions,
                                   idvar = "Sample",
                                   timevar = "Signature",
                                   direction = "wide")

# extract the relative contributions (the rescaled ones)
scaled_relative_contributions_df <- scaled_contributions_df[ ,seq(1, ncol(scaled_contributions_df), 2)]

# convert the counts of each trinucleotide to relative contributions per extracted signature
relative_contribution_trinucleotide_per_signature <- apply(fit_res_strict$signatures, 2, function(x) x/sum(x))

# extract the number of groups we have to perform the operation on
number_of_cluster <- nrow(scaled_contributions_df)

# create the empty list that will contain all the data.frames
relative_contribution_list <- c(rep(list(NA),
                                    number_of_cluster))

# add the names to the lists
names(relative_contribution_list) <- scaled_contributions_df[,1]

# for reach row (= sample/cluster) we create the relative contribution vector
for (i in 1:number_of_cluster){
  relative_contribution_vector <- as.numeric(scaled_contributions_df[i,seq(3, ncol(scaled_contributions_df), 2)]) 
  
  # mutliply vector * matrix
  mutliplied_matrix <- sweep(relative_contribution_trinucleotide_per_signature,
                             MARGIN=2,
                             relative_contribution_vector,
                             `*`)
  
  # store the matrix in the "master list"
  relative_contribution_list[[i]] <- mutliplied_matrix
}

# create the empty list that will contain all the data.frames
normalized_relative_contribution_list <- c(rep(list(NA),
                                               number_of_cluster))

# add the names to the lists
names(normalized_relative_contribution_list) <- scaled_contributions_df[,1]

# for reach row (= sample/cluster) we create the normalized relative contribution vector
for (i in 1:number_of_cluster){
  
  # get data frame
  normalized_relative_contribution_df <- relative_contribution_list[[i]] 
  
  # normalize the values of the data frame
  for (j in 1:nrow(normalized_relative_contribution_df)){
    vector_total <- sum(normalized_relative_contribution_df[j,])
    rescaled_vector <- as.numeric(normalized_relative_contribution_df[j, ])/vector_total
    normalized_relative_contribution_df[j, ] <- rescaled_vector
  }
  
  # save the normalized data frame
  normalized_relative_contribution_list[[i]] <- normalized_relative_contribution_df
}

# Function to reverse complement DNA sequence
reverse_complement_dna_sequence <- function(dnaSequence){
  return(reverse(chartr("ATGC", "TACG", dnaSequence)))
}

# Get probability for each trinucleotide context of selected SNVs
probability_df <- data.frame()
for (i in 1:nrow(SNV_list)){ 
  
  # Extract SNV info
  snv_line <- SNV_list[i,]
  snv <- paste0(snv_line$chr, ":", snv_line$start, "_", snv_line$ref, "/", snv_line$alt)
  
  # Get corresponding sample and cluster info
  snv_samplename <- SNV_list$sample_name[i]
  
  # Get normalized relative contributions belonging to the cluster
  if (snv_samplename %in% names(normalized_relative_contribution_list)){
    sample_cluster_probabilities <- normalized_relative_contribution_list[[which(names(normalized_relative_contribution_list) ==
                                                                                   snv_samplename)]]
  } else {
    sample_cluster_probabilities <- matrix(rep(rep(NA, 672)),nrow = 96,
                                           dimnames = list(rownames(normalized_relative_contribution_list[[1]]),
                                                           colnames(normalized_relative_contribution_list[[1]]))) %>%
      as.data.frame()
  }
  
  # Get mutation type and context for the mutation
  mutation_type <- paste0(snv_line$ref, ">", snv_line$alt)
  grl <- makeGRangesFromDataFrame(snv_line[,c("chr", "start", "end", "ref", "alt")],
                                  keep.extra.columns = T)
  GenomeInfoDb::genome(grl) = 'hg38'
  mut_context <- mut_context(grl, ref_genome = "hg38")
  
  # Get 96-nucleotide context of the mutation
  valid_mutations <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (mutation_type %in% valid_mutations){
    mut_context_96 <- paste0(substring(mut_context, first = 1, last = 1), 
                             "[", mutation_type, "]",
                             substring(mut_context, first = 3, last = 3))
    # Get reverse compliment if mutation type is on the negative strand
  } else {
    mut_context_96 <- paste0(reverse_complement_dna_sequence(substring(mut_context, first = 3, last = 3)),
                             "[", reverse_complement_dna_sequence(substring(mutation_type, first = 1, last = 1)), ">",
                             reverse_complement_dna_sequence(substring(mutation_type, first = 3, last = 3)), "]",
                             reverse_complement_dna_sequence(substring(mut_context, first = 1, last = 1)))
  }
  # Get mutation probabilities for the mutation context and add to data frame
  mut_context_probabilities <- sample_cluster_probabilities[mut_context_96,, drop = F]
  snv_line_compl <- cbind(snv_line, 
                          data.frame(trinucleotide_context = mut_context,
                                     trinucleotide_mutation = mut_context_96),
                          mut_context_probabilities)
  probability_df <- rbind(probability_df, snv_line_compl)
}

#Save probabilities
write_tsv(probability_df,
          "~/Projects/hypermutated_ALL/ANALYSES/pathogenicVariants/20251217_SNVprobabilities_pathogenic.tsv")
