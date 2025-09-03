#Load packages
library(tidyverse)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(MutationalPatterns)

# load the HG38 reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#Load vcf files
list_of_vcfs <- list.files("path/to/vcfs/",
                           pattern = ".vepAnnotated.vcf.gz$", full.names = T)

#Make empty matrix to which we will add the mutation matrices
mut_mat_collection_dbs <- NULL

#Loop over directories and make mutation matrix
for (vcf_file in list_of_vcfs){
  #Get patient id
  patient_id <- str_remove(basename(vcf_file), "_.+")
  #Load and filter vcf
  mnp_vcf <- loadVcfFile(vcf_file,
                         "BSgenome.Hsapiens.UCSC.hg38",
                         "chromosomes")
  GenomeInfoDb::genome(mnp_vcf)  <- "hg38"
  mnp_vcf_no_centromeric_variants <- excludeVariantsInCentromericRegions(mnp_vcf,
                                                                         "~/Projects/hypermutated_ALL/CODE/all_relapse_study/R/centromeres/cytoBand.hg38.centromeresOnly.txt")
  mnp_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(mnp_vcf_no_centromeric_variants,
                                                           0.01)
  mnp_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(mnp_vcf_gnomad_filtered,
                                                       0.01)
  tumor_sample_id_vcf_line <- grep("##tumor_sample=",
                                   readLines(vcf_file),
                                   value = TRUE)
  tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
  df_gonl_filtered <- convertVcfObjectToDataFrame(mnp_vcf_gonl_filtered)
  read_count_ref_idx <- which(names(df_gonl_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_ref"))
  read_count_alt_idx <- which(names(df_gonl_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_alt"))
  AF_idx <- which(names(df_gonl_filtered) %in% 
                    paste0(tumor_sample_id, "_read_based_AF"))
  df_gonl_filtered[paste0(tumor_sample_id, "_total_read_count")] <- 
    df_gonl_filtered[read_count_ref_idx] + 
    df_gonl_filtered[read_count_alt_idx]
  total_read_count_idx <- which(names(df_gonl_filtered) %in% 
                                  paste0(tumor_sample_id, "_total_read_count"))
  df_read_based_filtered <- df_gonl_filtered[df_gonl_filtered[,total_read_count_idx] >= 20 &
                                               df_gonl_filtered[,read_count_alt_idx] >= 5 &
                                               df_gonl_filtered[,AF_idx] >= 0.15,]
  mnp_filtered_vcf_object <- mnp_vcf[rownames(mnp_vcf) %in% 
                                       rownames(df_read_based_filtered), ]
  #Make mutation matrix
  sample_dbs_grl <- get_mut_type(rowRanges(mnp_filtered_vcf_object), "dbs", predefined_dbs_mbs=T)
  sample_dbs_grl <- get_dbs_context(sample_dbs_grl)
  sample_mut_mat <- count_dbs_contexts(sample_dbs_grl)
  colnames(sample_mut_mat) <- patient_id
  #Add mutation matrix to large matrix
  mut_mat_collection_dbs <- cbind(mut_mat_collection_dbs, sample_mut_mat)
}

#Save mutation matrix
subset_mut_mat <- mut_mat_collection_dbs[, colnames(mut_mat_collection_dbs) != ""]
mnp_mutation_matrix_ALCL = subset_mut_mat
save(object = mnp_mutation_matrix_ALCL,
     file = "/path/to/mnp_mutation_matrix.rdata")


