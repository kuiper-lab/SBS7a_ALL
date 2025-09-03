args = commandArgs(trailingOnly=TRUE)

#Load packages
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(MutationalPatterns)
library(ensemblVEP)
library(tidyverse)

# install the hg38 known gene info
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

# load the HG38 reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# load the hg38 genes
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

#Perform vcf filtering
sample_vcf_file <- args[1]
outputdir <- args[2]
sample_name <- gsub(".vcf.gz", "", basename(sample_vcf_file))

# test if the needed arguments are provided, if not, return an error
if (length(args)==0) {
  stop("Input file and output folder must be provided", call.=FALSE)
}

#Load vcf file
sample_vcf <- loadVcfFile(sample_vcf_file,
                          "BSgenome.Hsapiens.UCSC.hg38",
                          "chromosomes")


GenomeInfoDb::genome(sample_vcf)  <- "hg38"

#Remove mutations in centromeric regions
snp_vcf_no_centromeric_variants <- excludeVariantsInCentromericRegions(sample_vcf, "/hpc/pmc_kuiper/HypermutatedALL_project/REFERENCE/cytoBand.hg38.centromeresOnly.txt")

#Filter on population frequency
snp_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(snp_vcf_no_centromeric_variants,
                                                         0.01)

snp_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(snp_vcf_gnomad_filtered,
                                                     0.01)

#Apply read/vaf filter
tumor_sample_id_vcf_line <- grep("##tumor_sample=",
                                 readLines(sample_vcf_file),
                                 value = TRUE)
tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
df_gonl_filtered <- convertVcfObjectToDataFrame(snp_vcf_gonl_filtered)
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
snp_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                        rownames(df_read_based_filtered), ]
vcf_output <-  paste0(outputdir, "/", sample_name, "_filtered.vcf")
print("saving filtered vcf file to")
print(vcf_output)
writeVcf(snp_filtered_vcf_object, vcf_output)

#Calculate transcriptional strand bias
strand <- mut_strand(vcf = snp_filtered_vcf_object,
                     genes_hg38)
mut_mat_s <- mut_matrix_stranded(rowRanges(snp_filtered_vcf_object),
                                 ref_genome,
                                 genes_hg38)
colnames(mut_mat_s) <- tumor_sample_id
strand_counts <- strand_occurrences(mut_mat_s)
strand_counts$group <- tumor_sample_id

strand_bias <- strand_bias_test(strand_counts,
                                fdr_cutoffs = c(0.05, 0.01, 0.005),
                                p_cutoffs = c(0.05, 0.01, 0.005))


#Save files
print("Saving Transcriptional Strand Bias")
TSC_output <-  paste0(outputdir, "/", sample_name, "_strandcounts.csv")
write_csv(strand_counts, TSC_output)
TSB_output <- paste0(outputdir, "/", sample_name, "_strandbias.csv")
write_csv(strand_bias, TSB_output)

print("Ready")
