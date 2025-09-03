#Load packages
library(data.table)
library(scales)
library(cowplot)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(tidyverse)
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(readxl)

#Make function that identifies SBSs which overlap with CNAs
snvs_in_gained_region <- function(snvs, cnvs){
  cnvs$COPY_NUMBER <- cnvs$MEAN_COPY_RATIO * 2
  chr21_snvs <- snvs[snvs$chr == "chr21",]
  if ("chr21" %in% cnvs$CONTIG & nrow(chr21_snvs) > 0){
    chr21_cnvs <- cnvs[cnvs$CONTIG == "chr21",]
    for (variant in 1:nrow(chr21_snvs)){
      chromosome <- "chr21"
      position <- as.numeric(chr21_snvs$start[variant])
      type <- 'NO CNA'
      in_cnv <- FALSE
      mean_copy_ratio <- 1
      copy_number <- 2
      if (position >= min(chr21_cnvs$START)){
        closest_start <- which.max(chr21_cnvs$START[chr21_cnvs$START <= position])
        if (position >= chr21_cnvs$START[closest_start] &
            position <= chr21_cnvs$END[closest_start]){
          in_cnv <- TRUE
          mean_copy_ratio <- chr21_cnvs$MEAN_COPY_RATIO[closest_start]
          type <- chr21_cnvs$CALL[closest_start]
          copy_number <- chr21_cnvs$COPY_NUMBER[closest_start]
        }
      }
      chr21_snvs$IN_CNV[variant] <- in_cnv
      chr21_snvs$Type[variant] <- type
      chr21_snvs$MEAN_COPY_RATIO[variant] <- mean_copy_ratio
      chr21_snvs$COPY_NUMBER[variant] <- copy_number
    }
  } else if (nrow(chr21_snvs) != 0){
    chr21_snvs$IN_CNV <- FALSE
    chr21_snvs$Type <- 'NO CNA'
    chr21_snvs$MEAN_COPY_RATIO <- 1
    chr21_snvs$COPY_NUMBER <- 2
  } else {
    snvs$IN_CNV <- FALSE
    snvs$Type <- 'NO CNA'
    snvs$MEAN_COPY_RATIO <- 1
    snvs$COPY_NUMBER <- 2
    chr21_snvs <- snvs[snvs$chr == "chr21",]
  }
  return(chr21_snvs[chr21_snvs$Type == "+",])
}

#Get metadata
metadata <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/UV_germline_samples_overview.xlsx")

#Get list of vcfs
somatic_vcfs <- list.files("~/Projects/hypermutated_ALL/RESULTS/filteredVcf_multipletimepoints/VEP105_nogermlinefilter_nosoftclippedbases_010/",
                        full.names = T)

#Get list of cnv files
cnv_files <- list.files("~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/GATK_calls/",
                        full.names = T)
#Link files
rm(gained_somatic_snvs_combined)
for (patient in metadata$patient_id){
  #Get IDs
  tumor_id <- metadata$Tumor_ID[metadata$patient_id == patient]
  normal_id <- metadata$Normal_ID[metadata$patient_id == patient]
  subtype <- metadata$Subtype[metadata$patient_id == patient]
  
  #Read CNV file
  cnv_path <- grep(tumor_id, cnv_files, value = T)
  cnv_df <- read_tsv(cnv_path)
  cnv_df$MEAN_COPY_RATIO <- 2^cnv_df$MEAN_LOG2_COPY_RATIO
  
  #Read SNV file
  snv_vcf <- grep(paste0(patient, "Dx"), somatic_vcfs, value = T)
  snv_df <- convertVcfObjectToDataFrame(loadVcfFile(snv_vcf,
                                                    "BSgenome.Hsapiens.UCSC.hg38",
                                                    "chromosomes"))
  
  #Select SNVs in gains
  gained_snvs <- snvs_in_gained_region(snv_df, cnv_df)
  
  if(nrow(gained_snvs) >0){
    #Save in large dataframe
    if (sum(grepl("Dx", colnames(gained_snvs))) > 0){
      tumor_id <- paste0(patient, "Dx")
    }
    gained_snvs$depth <- gained_snvs[, paste0(tumor_id, "_read_count_alt")] + gained_snvs[, paste0(tumor_id, "_read_count_ref")]
    gained_snvs <- gained_snvs[, c(paste0(tumor_id, "_read_count_alt"),
                                   "depth",
                                   paste0(tumor_id, "_read_based_AF"),
                                   "chr",
                                   "start",
                                   "ref",
                                   "alt",
                                   "IN_CNV",
                                   "Type",
                                   "MEAN_COPY_RATIO",
                                   "COPY_NUMBER")]
    colnames(gained_snvs) <- c("alt_reads",
                               "depth",
                               "VAF",
                               "chr",
                               "start",
                               "ref",
                               "alt",
                               "IN_CNV",
                               "Type",
                               "MEAN_COPY_RATIO",
                               "COPY_NUMBER")
    gained_snvs$sample <- patient
    gained_snvs$subtype <- subtype
    if (!exists("gained_somatic_snvs_combined")){
      gained_somatic_snvs_combined <- gained_snvs
    } else {
      gained_somatic_snvs_combined <- rbind (gained_somatic_snvs_combined, gained_snvs)
    }
  }
}

#Save
write.table(gained_somatic_snvs_combined,
            "~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/SomaticGainedSNVs_UV_010.tsv",
            sep = "\t")