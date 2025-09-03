#Load packages
library(tidyverse)
library(MutationalPatterns)
library(VariantAnnotation)
library(stringi)

#Get list of vcfs
vcf_list <- list.files("~/Projects/SingleCellWGS/DATA/PTATO/Output_P0625/snvs/P0625/",
                       recursive = T,
                       pattern = ".snvs.ptato.filtered.vcf.gz$",
                       full.names = T)

#Go through vcfs and perform filtering
for (vcf_file in vcf_list){
  #Get sample ID
  sample <- basename(vcf_file) %>%
    str_remove(".snvs.ptato.filtered.vcf.gz") %>%
    str_remove("P0625.ann_")
  #Read in vcf
  vcf <- VariantAnnotation::readVcf(vcf_file)
  #Get vaf and alternative reads
  vcf_AD <- geno(vcf)$CAD %>% as.data.frame()
  vcf_alt <- apply(vcf_AD, 2, stri_extract_last, regex = "\\d+")
  vcf_alt_numeric <- apply(vcf_alt, 2, as.numeric)
  vcf_vaf <- geno(vcf)$VAF
  #Perform filtering
  vcf_filtered <- vcf[vcf_alt_numeric >= 5 &
                        vcf_vaf >= 0.25 &
                        info(vcf)$MQ >= 59]
  writeVcf(vcf_filtered,
           paste0("~/Projects/SingleCellWGS/ANALYSIS/FilterPTAvcf/P0625/PTATOvcf_adjustedVAF025andMQfiltered/",
                  sample,
                  ".vcf"))
}
