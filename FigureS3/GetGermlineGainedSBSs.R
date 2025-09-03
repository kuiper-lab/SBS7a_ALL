#Load packages
library(tidyverse)
library(MutationalPatterns)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(VariantAnnotation)
library(Polychrome)
library(readxl)

cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

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

#List data
metadata <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/UV_germline_samples_overview.xlsx")
haplotypecaller_list <- list.files("~/Projects/hypermutated_ALL/DATA/hg38/HaplotypeCaller/",
                                   pattern = "-HaplotypeCallerCalls.filter.passFiltered.snp.chr21.vcf.gz$",
                                   full.names = T,
                                   recursive = T)
cnvs <- list.files("~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/GATK_calls/",
                   full.names = T)

#Check availability of CNA data
for (tumor_id in metadata$Tumor_ID){
  print(tumor_id)
  print(grep(tumor_id, cnvs, value = T))
}

#Go through vcfs and filter
rm(gained_germline_snvs_combined)
for (patient in metadata$patient_id){
  tumor_id <- metadata$Tumor_ID[metadata$patient_id == patient]
  normal_id <- metadata$Normal_ID[metadata$patient_id == patient]
  subtype <- metadata$Subtype[metadata$patient_id == patient]
  
  #Load normal vcf
  normal_vcf_path <- grep(normal_id, haplotypecaller_list, value = T)
  normal_vcf <- loadVcfFile(normal_vcf_path,
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "chromosomes")
  GenomeInfoDb::genome(normal_vcf)  <- "hg38"
  
  #Load tumor vcf
  tumor_vcf_path <- grep(tumor_id, haplotypecaller_list, value = T)
  tumor_vcf <- loadVcfFile(tumor_vcf_path,
                           "BSgenome.Hsapiens.UCSC.hg38",
                           "chromosomes")
  GenomeInfoDb::genome(tumor_vcf)  <- "hg38"
  
  #Link to CNA data
  cnv_path <- grep(tumor_id, cnvs, value = T)
  
  #Filter out centromeric variants
  tumor_hg38_no_centromeric_variants <- excludeVariantsInCentromericRegions(tumor_vcf,
                                                                            "~/Projects/hypermutated_ALL/CODE/R/centromeres/cytoBand.hg38.centromeresOnly.txt")
  normal_hg38_no_centromeric_variants <- excludeVariantsInCentromericRegions(normal_vcf,
                                                                             "~/Projects/hypermutated_ALL/CODE/R/centromeres/cytoBand.hg38.centromeresOnly.txt")
  #Get overlapping variants
  tumor_hg38_no_centromeric_variants_overlap <- tumor_hg38_no_centromeric_variants[names(tumor_hg38_no_centromeric_variants) %in% names(normal_hg38_no_centromeric_variants)]
  normal_hg38_no_centromeric_variants_overlap <- normal_hg38_no_centromeric_variants[names(normal_hg38_no_centromeric_variants) %in% names(tumor_hg38_no_centromeric_variants)]
  
  #Filter on depth, vaf and alternative reads in normal
  AD <- geno(normal_hg38_no_centromeric_variants_overlap)$AD[,1]
  AD_df <- data.frame(unlist(lapply(AD, "[[", 1)), unlist(lapply(AD, "[[", 2)))
  colnames(AD_df) <- c("REF", "ALT")
  AD_df$DP <- geno(normal_hg38_no_centromeric_variants_overlap)$DP[,1]
  AD_df$VAF <- AD_df$ALT / AD_df$DP
  filtered_variants <- AD_df[AD_df$ALT >= 5 &
                               AD_df$DP >= 20 &
                               AD_df$VAF >= 0.25 & AD_df$VAF <= 0.75,]
  
  #Save filtered vcfs
  writeVcf(normal_hg38_no_centromeric_variants_overlap[names(normal_hg38_no_centromeric_variants_overlap) %in% rownames(filtered_variants)],
           paste0("~/Projects/hypermutated_ALL/DATA/hg38/HaplotypeCaller/UV_20082024_filtered/P_", patient, "F1.chr21.germlinefiltered.vcf"))
  writeVcf(tumor_hg38_no_centromeric_variants_overlap[names(tumor_hg38_no_centromeric_variants_overlap) %in% rownames(filtered_variants)],
           paste0("~/Projects/hypermutated_ALL/DATA/hg38/HaplotypeCaller/UV_20082024_filtered/P_", patient, "Dx.chr21.germlinefiltered.vcf"))
  
  #Make dataframe of filtered tumor data
  tumor_hg38_no_centromeric_variants_overlap_filtered <- tumor_hg38_no_centromeric_variants_overlap[names(tumor_hg38_no_centromeric_variants_overlap) %in% rownames(filtered_variants)]
  tumor_filtered_ALT <- unlist(lapply(geno(tumor_hg38_no_centromeric_variants_overlap_filtered)$AD[, 1], "[[", 2))
  tumor_filtered_DP <- geno(tumor_hg38_no_centromeric_variants_overlap_filtered)$DP[,1]
  tumor_filtered_df <- data.frame(tumor_filtered_ALT, tumor_filtered_DP, tumor_filtered_ALT / tumor_filtered_DP)
  colnames(tumor_filtered_df) <- c("alt_reads", "depth", "VAF")
  
  #Filter on tumor vaf
  tumor_filtered_df_Dxfiltered <- tumor_filtered_df[tumor_filtered_df$alt_reads >= 5 &
                                                      tumor_filtered_df$depth >= 20 &
                                                      tumor_filtered_df$VAF >= 0.1,]
  write.table(tumor_filtered_df_Dxfiltered,
              paste0("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/UV_filteredgermlineSNVs_010/P_",
                     patient,
                     "_germlineSNVsTumor.tsv"),
              sep = "\t")
  
  #Open snv dataframe and make new columns
  snv_df <- tumor_filtered_df_Dxfiltered
  snv_df$chr <- str_remove(rownames(snv_df), ":.+")
  snv_df$start <- as.numeric(str_remove(str_remove(rownames(snv_df), ".+:"), "_.+"))
  snv_df$ref <- str_remove(str_remove(rownames(snv_df), ".+_"), "/.+")
  snv_df$alt <- str_remove(str_remove(rownames(snv_df), ".+_"), ".+/")
  
  #Open cnv dataframe
  cnv_df <- read_tsv(cnv_path)
  cnv_df$MEAN_COPY_RATIO <- 2^cnv_df$MEAN_LOG2_COPY_RATIO
  
  #Select snvs in gains
  gained_snvs <- snvs_in_gained_region(snv_df, cnv_df)
  
  #Save in large dataframe
  gained_snvs$sample <- patient
  gained_snvs$subtype <- subtype
  if (!exists("gained_germline_snvs_combined")){
    gained_germline_snvs_combined <- gained_snvs
  } else {
    gained_germline_snvs_combined <- rbind(gained_germline_snvs_combined, gained_snvs)
  }
}

#Save
write.table(gained_germline_snvs_combined,
            "~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/GermlineGainedSNVs_UV_010.tsv",
            sep = "\t")