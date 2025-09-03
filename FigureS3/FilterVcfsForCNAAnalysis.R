#Load packages
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(MutationalPatterns)
library(tidyverse)
library(grid)
library(gridExtra)
library(stringi)
library(cowplot)
library(readxl)
library(cowplot)
library(grid)

# load the HG38 reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# install the hg38 known gene info
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

# load the hg38 genes
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

#Load cosmic signatures
cancer_signatures <- get_known_signatures()

#Make list of files
vcf_vector <- list.files("~/Projects/hypermutated_ALL/DATA/hg38/UV_WGS_nosoftclippedbases//",
                         pattern = "-merged-mutect2Calls.passFiltered.snp-WGS.vepAnnotated.vcf.gz$",
                         full.names = T)

#Filter the vcf files
for (vcf_path in vcf_vector) {
  #Extract patient id
  patient_id <- str_extract(vcf_path, "\\d\\d\\d\\d+")
  #Read vcf file
  sample_vcf <- loadVcfFile(vcf_path,
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  #Select patients with more than 2 samples
  if (length(colnames(sample_vcf)) > 2) {
    #Rename sample names
    for (sample_id in colnames(sample_vcf)) {
      if (endsWith(sample_id, patient_id)) {
        time <- str_split(sample_id, "-")[[1]][1]
        new_name <- paste0(patient_id, time)
      } else if (startsWith(sample_id, "P_")) {
        time <- str_split(sample_id, "_")[[1]][3]
        new_name <- paste0(patient_id, time)
      } else {
        new_name <- sample_id
      }
      if (endsWith(new_name, "D")) {
        new_name <- paste0(new_name, "x")
      }
      colnames(sample_vcf)[colnames(sample_vcf) == sample_id] <-
        new_name
    }
    #Filter out centromeric variants
    sample_vcf_no_centromeric_variants <-
      excludeVariantsInCentromericRegions(
        sample_vcf,
        "~/Projects/hypermutated_ALL/CODE/all_relapse_study/R/centromeres/cytoBand.hg38.centromeresOnly.txt"
      )
    #Filter on population frequency
    sample_vcf_gnomad_filtered <-
      filterOnGnomadAlleleFrequency(sample_vcf_no_centromeric_variants,
                                    0.01)
    sample_vcf_gonl_filtered <-
      filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                  0.01)
    sample_gonl_filtered_df <-
      convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
    
    #Make df with only tumor data
    timepoints <-
      grep("D|([^C]R)", colnames(sample_vcf_gonl_filtered), value = T)
    for (timepoint in 1:length(timepoints)) {
      time_name <- timepoints[timepoint]
      time_df <-
        sample_gonl_filtered_df[, grep(time_name, colnames(sample_gonl_filtered_df))]
      time_df <- time_df[, !grepl("mutect", colnames(time_df))]
      time_df[, paste0(time_name, "_read_count_total")] <-
        time_df[, paste0(time_name, "_read_count_ref")] +
        time_df[, paste0(time_name, "_read_count_alt")]
      if (timepoint == 1) {
        tumor_df <- time_df
      } else{
        tumor_df <- cbind(tumor_df, time_df)
      }
    }
    
    #Apply filter
    alt_cols <- grep("read_count_alt", colnames(tumor_df))
    total_cols <- grep("read_count_total", colnames(tumor_df))
    vaf_cols <- grep("read_based_AF", colnames(tumor_df))
    
    sample_filtered_df <-
      sample_gonl_filtered_df[(apply(tumor_df[, alt_cols], 1, max) >= 5) &
                                (apply(tumor_df[, total_cols], 1, max) >= 20) &
                                (apply(tumor_df[, vaf_cols], 1, max) >= 0.1), ]
    
    #Look at separate timepoints
    timepoints <-
      grep("D|([^C]R)", colnames(sample_vcf_gonl_filtered), value = T)
    for (timepoint in timepoints){
      read_count_ref_idx <- which(names(sample_filtered_df) %in% 
                                    paste0(timepoint, "_read_count_ref"))
      read_count_alt_idx <- which(names(sample_filtered_df) %in% 
                                    paste0(timepoint, "_read_count_alt"))
      AF_idx <- which(names(sample_filtered_df) %in% 
                        paste0(timepoint, "_read_based_AF"))
      sample_filtered_df[paste0(timepoint, "_total_read_count")] <- 
        sample_filtered_df[read_count_ref_idx] + 
        sample_filtered_df[read_count_alt_idx]
      total_read_count_idx <- which(names(sample_filtered_df) %in% 
                                      paste0(timepoint, "_total_read_count"))
      df_read_based_filtered <- sample_filtered_df[sample_filtered_df[,total_read_count_idx] >= 20 &
                                                     sample_filtered_df[,read_count_alt_idx] >= 5 &
                                                     sample_filtered_df[,AF_idx] >= 0.1,]
      sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                                 rownames(df_read_based_filtered), ]
      writeVcf(
        sample_filtered_vcf_object,
        paste0(
          "~/Projects/hypermutated_ALL/RESULTS/filteredVcf_multipletimepoints/VEP105_nogermlinefilter_nosoftclippedbases_010/P_",
          timepoint,
          ".vcf"
        )
      )
    }
  #Select patients with 2 timepoints
  } else if (length(colnames(sample_vcf)) == 2){
    #Filter out centromeric variants
    sample_vcf_no_centromeric_variants <-
      excludeVariantsInCentromericRegions(
        sample_vcf,
        "~/Projects/hypermutated_ALL/CODE/all_relapse_study/R/centromeres/cytoBand.hg38.centromeresOnly.txt"
      )
    #Filter on population frequency
    sample_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(sample_vcf_no_centromeric_variants,
                                                                0.01)
    sample_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                                            0.01)
    #Apply read/vaf filter
    df_gonl_filtered <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
    tumor_sample_id_vcf_line <- grep("##tumor_sample=",
                                     readLines(vcf_path),
                                     value = TRUE)
    tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
    if (grepl(patient_id, tumor_sample_id)){
      file_name <- tumor_sample_id
    } else {
      file_name <- paste0(patient_id, "Dx")
    }
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
                                                 df_gonl_filtered[,AF_idx] >= 0.1,]
    sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                               rownames(df_read_based_filtered), ]
    writeVcf(
      sample_filtered_vcf_object,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/filteredVcf_multipletimepoints/VEP105_nogermlinefilter_nosoftclippedbases_010/P_",
        file_name,
        ".vcf"
      )
    )
  }
}
