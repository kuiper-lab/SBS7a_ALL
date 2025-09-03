#Load packages
library(tidyverse)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(readxl)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse

#Get list of NER genes
NER_genes <- scan("~/Projects/hypermutated_ALL/ANALYSES/NER/KEGG_NER_symbols.txt",
                  character(), quote = "", sep = "\n")

#Get list of diagnosis vcfs
list_of_vcfs <- list.files("~/Projects/hypermutated_ALL/RESULTS/filteredVcf_multipletimepoints/VEP105_nogermlinefilter_nosoftclippedbases_015/",
                           pattern = "Dx", full.names = T)

#Go through vcfs
for (vcf_path in list_of_vcfs){
  #Read vcf file
  sample_vcf <- loadVcfFile(vcf_path,
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  #Make dataframe
  sample_df <- convertVcfObjectToDataFrame(sample_vcf)
  #Check if any mutations are left
  #Add sample name
  sample_df$sample <- basename(vcf_path) %>% str_remove(".vcf")
  #Add VAF
  AF_cols <- grep("read_based_AF", colnames(sample_df), value = T)
  if (length(AF_cols) > 2){
    diagnosis_AF <- grep("Dx", AF_cols, value = T)
  } else {
    diagnosis_AF <- names(which.max(colMeans(sample_df[, AF_cols])))
  }
  sample_df$VAF <- sample_df[, diagnosis_AF]
  #Add to larger dataframe
  sample_df_minimal <- sample_df[, c("chr", "start", "end", "ref", "alt",
                                     "Consequence", "IMPACT", "SYMBOL",
                                     "CADD_PHRED", "sample", "VAF")]
  if (vcf_path == list_of_vcfs[1]){
    combined_df <- sample_df_minimal
  } else {
    combined_df <- rbind(combined_df, sample_df_minimal)
  }
}

#Check if there are any mutations in NER
combined_df_NER <- combined_df[combined_df$SYMBOL %in% NER_genes,]

#Save
write.table(combined_df_NER,
            "~/Projects/hypermutated_ALL/ANALYSES/NER/AllNERMutations.tsv",
            sep = "\t")
