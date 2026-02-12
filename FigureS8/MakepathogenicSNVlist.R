#Load packages
library(tidyverse)
source("~/git/pmc_kuiper_projects/HypermutatedALL/CODE/local/createdFunctions_HM-ALL.R")
library(readxl)

#Load metadata
metadata_SBS7a <- read_excel("~/Projects/hypermutated_ALL/DOCS/SBS7aCohort_metadata.xlsx")
driver_gene_table <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/pathogenicVariants/41375_2024_2403_MOESM2_ESM.xlsx",
                                sheet = 1, skip = 3)
driver_list <- driver_gene_table$Gene

#Get list of vcfs
vcf_list <- list.files("~/Projects/hypermutated_ALL/RESULTS/filteredClusteredVcf/VEP105_nogermlinefilter_nosoftclippedbases/",
                       full.names = T)

#Get list of cluster files
cluster_list <- list.files( "~/Projects/hypermutated_ALL/ANALYSES/SBS7a/clustered_SNVs_nogermlinefilter_nosoftclippedbases_23042025/",
                            full.names = T)

#Get list of anonymous IDs
metadata <- as.data.frame(read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx"))
recode_vector <- setNames(metadata$Coded_num, metadata$SKION)
recode_vector <- recode_vector[recode_vector %in% metadata_SBS7a$`Anonymous ID`[!is.na(metadata_SBS7a$`Anonymous ID`)]]

#Make list of all SNVs
for (snv_vcf in vcf_list){
  #Get patient ID
  patient_id <- snv_vcf %>% basename() %>% str_remove(".vcf")
  #Load vcf file
  sample_vcf <- loadVcfFile(snv_vcf,
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  #Get some information
  sample_df <- convertVcfObjectToDataFrame(sample_vcf)
  minimal_df <- data.frame(patient = patient_id,
                           chr = sample_df$chr,
                           start = sample_df$start,
                           end = sample_df$end,
                           ref = sample_df$ref,
                           alt = sample_df$alt,
                           CADD = as.numeric(sample_df$CADD_PHRED),
                           gene = sample_df$SYMBOL,
                           protein_position = sample_df$Protein_position,
                           amino_acid = sample_df$Amino_acids,
                           anonymous = str_replace_all(str_extract(patient_id,
                                                                   "\\d\\d\\d\\d+"),
                                                       recode_vector) %>%
                             str_remove("P_"))
  
  #Get cluster information
  cluster_file <- grep(patient_id,
                       cluster_list,
                       value = T)
  if (length(cluster_file) > 0){
    cluster_info <- read.table(cluster_file)
    minimal_df <- minimal_df[names(sample_vcf) %in% rownames(cluster_info),]
    cluster_info_sorted <- cluster_info[names(sample_vcf)[names(sample_vcf) %in% rownames(cluster_info)],]
    minimal_df$cluster <- cluster_info_sorted$cluster
    minimal_df$sample_name <- paste0(minimal_df$anonymous, "_cluster", minimal_df$cluster)
  } else {
    minimal_df$cluster <- "Dx"
    minimal_df$sample_name <- paste0(minimal_df$anonymous, "Dx")
  }
  if (snv_vcf == vcf_list[1]){
    combined_df <- minimal_df
  } else {
    combined_df <- rbind(combined_df,
                         minimal_df)
  }
}

#Select pathogenic mutations in driver genes
combined_df_pathogenic <- combined_df[combined_df$CADD >= 15,]
combined_df_driver <- combined_df_pathogenic[combined_df_pathogenic$gene %in% driver_list,]

#Save SNV table
write.table(combined_df_driver,
            "~/Projects/hypermutated_ALL/ANALYSES/pathogenicVariants/20251217_clustered_SNV_list_pathogenic.tsv",
            sep = "\t")

