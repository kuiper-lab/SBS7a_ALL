#Load packages
library(MutationalPatterns)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(tidyverse)
library(grid)
library(gridExtra)
library(stringi)
library(cowplot)
library(readxl)
library(cowplot)
library(grid)
library(plyr)
library(lsa)

PLOT_AF_AND_PROFILE_PLOT=TRUE
WRITE_VARIANTS_TO_FILE=T

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
vcf_vector <- list.files("~/Projects/hypermutated_ALL/DATA/hg38/UV_WGS_nosoftclippedbases/",
                         pattern = "-merged-mutect2Calls.passFiltered.snp-WGS.vepAnnotated.vcf.gz$",
                         full.names = T)

#Make empty matrix to which we will add the mutation matrices
mut_mat_collection <- matrix(nrow=96)

mut_mat_collection <- NULL

#Filter the vcf files and make mutation matrix for each timepoint in each patient
rm(combined_strand_counts)
rm(combined_strand_bias)
for (vcf_path in vcf_vector) {
  #Get patient id
  patient_id <- str_extract(vcf_path, "\\d\\d\\d\\d+")
  #Load vcf
  sample_vcf <- loadVcfFile(vcf_path,
                            "BSgenome.Hsapiens.UCSC.hg38",
                            "chromosomes")
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  #Take patients with more than 2 samples
  if (length(colnames(sample_vcf)) > 2) {
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
      sample_gonl_filtered_df[(apply(tumor_df[, alt_cols], 1, max) >= 3) &
                                (apply(tumor_df[, total_cols], 1, max) >= 20) &
                                (apply(tumor_df[, vaf_cols], 1, max) >= 0.25), ]
    sample_hg38_filtered <-
      sample_vcf_gonl_filtered[rownames(sample_filtered_df)]
    #Save vcf
    writeVcf(
      sample_hg38_filtered,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/filteredClusteredVcf/VEP105_nogermlinefilter_nosoftclippedbases/P_",
        patient_id,
        ".vcf"
      )
    )
    
    #Get allele frequencies matrix
    sample_af_matrix <-
      extractAllelFrequenciesFromVcf(sample_hg38_filtered)
    
    #Select only read bases allele frequencies
    sample_af_df <-
      sample_af_matrix[, (seq(ncol(sample_af_matrix) / 4) * 4) - 1]
    loginfo(paste0(
      "Using the read based allele frequencies of the vcf-file to cluster variants"
    ))
    
    #Select tumors
    sample_tumour_ids <- grep("D|([^C]R)", colnames(sample_af_df))
    sample_af_df_tumours <- sample_af_df[, sample_tumour_ids]
    sample_transformed_tumours <- sample_af_df_tumours
    
    #Transform allele frequencies to binary
    for (colname in colnames(sample_transformed_tumours)) {
      sample_transformed_tumours[, colname] = ifelse(sample_af_df_tumours[, colname] > 0, 1, 0)
      
    }
    
    #Get the different patterns in the data
    sample_countdf <- plyr::count(sample_transformed_tumours)
    sample_significant_patterns <- sample_countdf
    
    #Cluster mutations
    clusters <- c()
    for (AF_row in 1:nrow(sample_transformed_tumours)) {
      cluster <- NA
      cossim_score <- 0
      for (pattern_row in 1:nrow(sample_significant_patterns)) {
        if (identical(
          as.numeric(sample_transformed_tumours[AF_row,]),
          as.numeric(sample_significant_patterns[pattern_row, 1:ncol(sample_significant_patterns) -
                                                 1])
        )) {
          cluster <- pattern_row
          clusters <- append(clusters, cluster)
          break
        }
      }
    }
    sample_transformed_tumours_clustered <-
      as.data.frame(sample_af_df_tumours)
    sample_transformed_tumours_clustered$cluster <- clusters
    sample_filtered_variants_clusters_df <-
      sample_transformed_tumours_clustered
    
    #Count mutations per sample
    for (colname in colnames(sample_filtered_variants_clusters_df)) {
      if (colname != "cluster") {
        counts <- sum(sample_filtered_variants_clusters_df[, colname] > 0)
        print(paste(colname, counts, sep = ","))
      }
    }
    
    #Rename clusters
    unique(sample_filtered_variants_clusters_df$cluster)
    
    #Save clustered dataframe
    write.csv(
      sample_filtered_variants_clusters_df,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/ClusteredVariantsVafs_nogermlinefilter_nosoftclippedbases/P_",
        patient_id,
        ".csv"
      )
    )
    save(
      object = sample_filtered_variants_clusters_df,
      file = paste0(
        "~/Projects/hypermutated_ALL/RESULTS/ClusteredVariantsVafs_nogermlinefilter_nosoftclippedbases/P_",
        patient_id,
        ".rdata"
      )
    )
    #Split where necessary
    if (patient_id == {PATIENT}) {
      sample_filtered_variants_clusters_df <- splitClusterOnTimePointSpecificAlleleFrequency(cluster_df = sample_filtered_variants_clusters_df,
                                                                                             time_point = {TIME_POINT},
                                                                                             cluster = {CLUSTER},
                                                                                             AF_threshold = {VAF_THRESHOLD})
    }
    
    #Count number of mutations per cluster
    cluster_count <-
      plyr::count(sample_filtered_variants_clusters_df$cluster)
    
    #Remove small clusters (less than 75 mutations)
    sample_small_clusters <-
      cluster_count[!(cluster_count$freq > 75), ]
    sample_filtered_variants_clusters_df <-
      removeClustersFromDataFrame(cluster_data_frame = sample_filtered_variants_clusters_df,
                                  clusters_to_remove = sample_small_clusters$x)
    
    sample_clusterSpecificAlleleFrequencyPlots <-
      plotClusterSpecificAlleleFrequency(alleleFrequencyDataFrame = sample_filtered_variants_clusters_df,
                                         minimal = TRUE)
    
    #split the vcf based on the clustering
    sample_filtered_variants_clusters <-
      splitVcfBasedOnCluster(sample_hg38_filtered,
                             sample_filtered_variants_clusters_df)
    
    #Add usefull headers
    names(sample_filtered_variants_clusters) <-
      paste0(
        patient_id,
        "_" ,
        names(sample_filtered_variants_clusters)
      )
    
    # create a mutation matrix
    sample_mut_mat <-
      mut_matrix(sample_filtered_variants_clusters , ref_genome = ref_genome)
    write.csv(
      sample_mut_mat,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/ClusteredMutationMatrices_nogermlinefilter_nosoftclippedbases_23042025/P_",
        patient_id,
        ".csv"
      )
    )
    #Split vcf based on clusters
    split_vcf <- splitVcfObjectOnClusters(sample_hg38_filtered,
                                          sample_filtered_variants_clusters_df)
    names(split_vcf) <- paste0(
      patient_id,
      "_cluster" ,
      names(split_vcf)
    )
    #Calculate transcriptional strand bias
    for (vcf_name in names(split_vcf)){
      strand <- mut_strand(vcf = split_vcf[[vcf_name]],
                           genes_hg38)
      mut_mat_s <- mut_matrix_stranded(rowRanges(split_vcf[[vcf_name]]),
                                       ref_genome,
                                       genes_hg38)
      colnames(mut_mat_s) <- vcf_name
      strand_counts <- strand_occurrences(mut_mat_s)
      strand_counts$group <- vcf_name
      if (exists("combined_strand_counts")){
        combined_strand_counts <- rbind(combined_strand_counts, strand_counts)
      } else {
        combined_strand_counts <- strand_counts
      }
      strand_bias <- strand_bias_test(strand_counts,
                                      fdr_cutoffs = c(0.05, 0.01, 0.005),
                                      p_cutoffs = c(0.05, 0.01, 0.005))
      
      if (exists("combined_strand_bias")){
        combined_strand_bias <- rbind(combined_strand_bias, strand_bias)
      } else {
        combined_strand_bias <- strand_bias
      }
    }
    #Plot clusters
    if (PLOT_AF_AND_PROFILE_PLOT) {
      sample_af_merged_plots <-
        plot_grid(plotlist = sample_clusterSpecificAlleleFrequencyPlots,
                  ncol = 1)
      
      sample_profile_plots <- plot_96_profile(sample_mut_mat,
                                              ymax = 0.30,
                                              condensed = T)
      
      fn <-
        paste0(
          "~/Projects/hypermutated_ALL/RESULTS/ClusterVafPlots_nogermlinefilter_nosoftclippedbases/P_",
          patient_id,
          "-",
          "overview",
          ".pdf"
        )
      pdf(file = fn)
      multiplot(sample_af_merged_plots,
                sample_profile_plots,
                layout = matrix(c(1, 2), ncol = 2))
      dev.off()
      
      cancer_signatures_sbs2_sbs13 <-
        cbind(cancer_signatures, rowSums(cancer_signatures[, c("SBS2", "SBS13")]))
      colnames(cancer_signatures_sbs2_sbs13)[length(colnames(cancer_signatures_sbs2_sbs13))] <-
        "SBS2/SBS13"
      cosmic_order_sbs2_sbs13 = colnames(cancer_signatures_sbs2_sbs13)
      cos_sim_samples_signatures = cos_sim_matrix(sample_mut_mat,
                                                  cancer_signatures_sbs2_sbs13)
      cosine_heatmap <-
        plot_cosine_heatmap(
          cos_sim_samples_signatures,
          cluster_rows = FALSE,
          plot_values = FALSE
        )
      
      ggsave(
        plot = cosine_heatmap,
        width = 10,
        height = 7,
        filename = paste0(
          "~/Projects/hypermutated_ALL/ANALYSES/MultipleRelapses/ClusterVafPlots/heatmaps/P_",
          patient_id,
          "-",
          "cosine_heatmap",
          ".jpg"
        )
      )
      
    }
    
    
    mut_mat_collection <- cbind(mut_mat_collection, sample_mut_mat)
    
    if (WRITE_VARIANTS_TO_FILE) {
      VARIANT_FILE_PATH <-
        paste0(
          "~/Projects/hypermutated_ALL/ANALYSES/SBS7a/clustered_SNVs_nogermlinefilter_nosoftclippedbases_23042025/P_",
          patient_id,
          "-clustered_SNVs.txt"
        )
      
      # write variants to tabular format:
      writeDataFrameToTablularFile(dataFrame = sample_filtered_variants_clusters_df,
                                   outputFile = VARIANT_FILE_PATH)
      
    }
    #Select patients with only 2 samples
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
    df_gonl_filtered <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
    #Filter on reads in tumor
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
                                                 df_gonl_filtered[,AF_idx] >= 0.15,]
    sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in% 
                                               rownames(df_read_based_filtered), ]
    #Save filtered vcf
    writeVcf(
      sample_filtered_vcf_object,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/filteredClusteredVcf/VEP105_nogermlinefilter_nosoftclippedbases/P_",
        patient_id,
        ".vcf"
      )
    )
    
    #Make mutation matrix
    sample_mut_mat <- mut_matrix(list(rowRanges(sample_filtered_vcf_object)),
                                 ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
    colnames(sample_mut_mat) <- file_name
    
    #Calculate transcriptional strand bias
    strand <- mut_strand(vcf = sample_filtered_vcf_object,
                         genes_hg38)
    mut_mat_s <- mut_matrix_stranded(rowRanges(sample_filtered_vcf_object),
                                     ref_genome,
                                     genes_hg38)
    colnames(mut_mat_s) <- file_name
    strand_counts <- strand_occurrences(mut_mat_s)
    strand_counts$group <- file_name
    if (exists("combined_strand_counts")){
      combined_strand_counts <- rbind(combined_strand_counts, strand_counts)
    } else {
      combined_strand_counts <- strand_counts
    }
    strand_bias <- strand_bias_test(strand_counts,
                                    fdr_cutoffs = c(0.05, 0.01, 0.005),
                                    p_cutoffs = c(0.05, 0.01, 0.005))
    
    if (exists("combined_strand_bias")){
      combined_strand_bias <- rbind(combined_strand_bias, strand_bias)
    } else {
      combined_strand_bias <- strand_bias
    }
    #Add sample mutation matrix to large mutatation matrix
    mut_mat_collection <- cbind(mut_mat_collection, sample_mut_mat)
    #Save mutation matrix
    write.csv(
      sample_mut_mat,
      paste0(
        "~/Projects/hypermutated_ALL/RESULTS/ClusteredMutationMatrices_nogermlinefilter_nosoftclippedbases_23042025/P_",
        patient_id,
        ".csv"
      )
    )
  }
}

#Save files
mut_mat_collection <- mut_mat_collection[, colnames(mut_mat_collection) != ""]
mut_mat_collection_UV_23042025 <- mut_mat_collection
save(object = mut_mat_collection_UV_23042025,
     file = "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/mut_mat_collection_UV_clustered_nogermlinefilter_nosoftclippedbases_23042025.rdata")
write.csv(mut_mat_collection_UV_23042025, "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/mut_mat_collection_UV_clustered_nogermlinefilter_nosoftclippedbases_23042025.csv")

save(object = combined_strand_counts,
     file = "~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/ALLstrandcounts_clustered_nogermlinefilter_nosoftclippedbases_23042025.rdata")
save(object = combined_strand_bias,
     file = "~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/ALLstrandbias_clustered_nogermlinefilter_nosoftclippedbases_23042025.rdata")


