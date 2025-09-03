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
library(readxl)
library(multimode)
library(diptest)
library(LaplacesDemon)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# load the HG38 reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#Load metadata
metadata <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/metadata.xlsx")
anonymous_ids <- read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx")

#Make list of files
vcf_vector <- list.files("~/Projects/hypermutated_ALL/DATA/hg38/UV_WGS_nosoftclippedbases//",
                         pattern = "-merged-mutect2Calls.passFiltered.snp-WGS.vepAnnotated.vcf.gz$",
                         full.names = T)
vcf_vector <- grep(paste0(metadata$Patient, collapse = "|"),
                   vcf_vector,
                   value = T)

#Make empty matrix to which we will add the mutation matrices
mut_mat_collection_combined <- NULL
vaf_df <- data.frame("Sample" = NA,
                     "VAF" = NA,
                     "Position" = NA,
                     "Modality" = NA,
                     "cutoff" = NA,
                     "subclonal_mutationload" = NA,
                     "clonal_mutationload" = NA)

#Filter the vcf files
for (vcf_path in vcf_vector) {
  #Extract skion number
  skion_number <- str_extract(vcf_path, "\\d\\d\\d\\d+")
  #Get anonymous ID
  anonymous_id <- anonymous_ids$Coded_num[anonymous_ids$SKION == skion_number]
  #Read vcf file
  sample_vcf_complete <- loadVcfFile(vcf_path,
                                     "BSgenome.Hsapiens.UCSC.hg38",
                                     "chromosomes")
  GenomeInfoDb::genome(sample_vcf_complete)  <- "hg38"
  #Get only chromosomes without CNVs
  CNV_blacklist <- metadata$CNVs[metadata$Patient == skion_number] %>%
    str_split(pattern = ",") %>%
    unlist()
  vcf_chromosomes <- names(sample_vcf_complete) %>%
    str_remove(":.+")
  sample_vcf <- sample_vcf_complete[!vcf_chromosomes %in% CNV_blacklist]
  
  #Rename sample names for samples with >2 timepoints
  if (length(colnames(sample_vcf)) > 2) {
    for (sample_id in colnames(sample_vcf)) {
      if (endsWith(sample_id, skion_number)) {
        time <- str_split(sample_id, "-")[[1]][1]
        new_name <- paste0(skion_number, time)
      } else if (startsWith(sample_id, "P_")) {
        time <- str_split(sample_id, "_")[[1]][3]
        new_name <- paste0(skion_number, time)
      } else {
        new_name <- sample_id
      }
      if (endsWith(new_name, "D")) {
        new_name <- paste0(new_name, "x")
      }
      colnames(sample_vcf)[colnames(sample_vcf) == sample_id] <-
        new_name
    }
  }
  #Rename sample names for samples with 2 timepoints
  if (length(colnames(sample_vcf)) == 2) {
    tumor_sample_id_vcf_line <- grep("##tumor_sample=",
                                     readLines(vcf_path),
                                     value = TRUE)
    tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
    colnames(sample_vcf)[colnames(sample_vcf) == tumor_sample_id] <- paste0(skion_number, "Dx")
    colnames(sample_vcf)[colnames(sample_vcf) != paste0(skion_number, "Dx")] <-  paste0(skion_number, "F1")
  }
  #We will now only look at Dx
  dx_vcf <- sample_vcf[, paste0(skion_number, "Dx")]
  
  #Filter out centromeric variants
  sample_vcf_no_centromeric_variants <-
    excludeVariantsInCentromericRegions(
      dx_vcf,
      "~/Projects/hypermutated_ALL/CODE/all_relapse_study/R/centromeres/cytoBand.hg38.centromeresOnly.txt"
    )
  #Filter on population frequency
  sample_vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(sample_vcf_no_centromeric_variants,
                                                              0.01)
  sample_vcf_gonl_filtered <- filterOnGonlAlleleFrequency(sample_vcf_gnomad_filtered,
                                                          0.01)
  #Apply read/vaf filter
  df_gonl_filtered <- convertVcfObjectToDataFrame(sample_vcf_gonl_filtered)
  read_count_ref_idx <- grep("_read_count_ref", names(df_gonl_filtered))
  read_count_alt_idx <- grep("_read_count_alt", names(df_gonl_filtered))
  AF_idx <- grep("_read_based_AF", names(df_gonl_filtered))
  df_gonl_filtered[paste0(colnames(dx_vcf), "_total_read_count")] <- 
    df_gonl_filtered[read_count_ref_idx] + 
    df_gonl_filtered[read_count_alt_idx]
  total_read_count_idx <- grep("_total_read_count", names(df_gonl_filtered))
  #Get all mutations and collect vaf distibutions
  df_read_based_filtered_all <- df_gonl_filtered[df_gonl_filtered[,total_read_count_idx] >= 20 &
                                                   df_gonl_filtered[,read_count_alt_idx] >= 5 &
                                                   df_gonl_filtered[,AF_idx] >= 0.15,]
  vaf_df_sample <- data.frame("Sample" = str_replace(colnames(dx_vcf),
                                                     skion_number,
                                                     anonymous_id),
                              "VAF" = df_read_based_filtered_all[,AF_idx],
                              "Position" = rownames(df_read_based_filtered_all))
  #Get VAF peak
  VAF_modes <- Modes(vaf_df_sample$VAF)$modes
  VAF_modes <- sort(VAF_modes, decreasing = T)
  vaf_peak <- max(VAF_modes)
  #Go to next if the VAF peak is below 0.25
  if (vaf_peak < 0.25){
    next()
  }
  #Test for modality
  if (is.unimodal(vaf_df_sample$VAF)){
    vaf_df_sample$Modality <- "unimodal"
  } else if (is.bimodal(vaf_df_sample$VAF)){
    vaf_df_sample$Modality <- "bimodal"
  } else if (is.multimodal(vaf_df_sample$VAF)){
    vaf_df_sample$Modality <- "multimodal"
  }
  #Check if there are multiple modes
  if (length(VAF_modes) > 1){
    for (mode_idx in 1:(length(VAF_modes)-1)){
      #Get modes
      high_mode <- max(VAF_modes[c(mode_idx, mode_idx + 1)])
      low_mode <- min(VAF_modes[c(mode_idx, mode_idx + 1)])
      #Get density
      sample_density <- density(vaf_df_sample$VAF)
      #Get density between modes
      modal_density_x <- sample_density$x[sample_density$x > low_mode &
                                            sample_density$x < high_mode]
      modal_density_y <- sample_density$y[sample_density$x > low_mode &
                                            sample_density$x < high_mode]
      vaf_antimode <- modal_density_x[which.min(modal_density_y)]
      #Get cutoffs and number of clonal mutations
      if (mode_idx == 1){
        cutoff <- vaf_antimode
        #Get number of clonal mutations
        clonal_mutationload <- sum(vaf_df_sample$VAF > vaf_antimode)
        #Add to vaf_df
        vaf_df_sample$clonal_mutationload <- clonal_mutationload
      } else {
        cutoff <- c(cutoff, vaf_antimode)
      }
    }
    vaf_df_sample$cutoff <- rep(list(cutoff), nrow(vaf_df_sample))
    #Get mutations from clone and all subclones
    cutoff_vector <- sort(unique(unlist(vaf_df_sample$cutoff)), decreasing = T)
    for (cutoff_idx in 1:(length(cutoff_vector) + 1)){
      #Filter mutations using the cutoff(s)
      if (cutoff_idx == 1){
        fraction <- "clonal"
        df_read_based_fraction <- df_read_based_filtered_all[vaf_df_sample$VAF > cutoff_vector[1],]
      } else {
        df_read_based_fraction <- df_read_based_filtered_all[vaf_df_sample$VAF < cutoff_vector[cutoff_idx-1] &
                                                               vaf_df_sample$VAF >= c(cutoff_vector, 0)[cutoff_idx],]
        fraction <- paste0("subclonal_", cutoff_idx - 1)
        #Get number of subclonal mutations
        subclonal_mutationload <- nrow(df_read_based_fraction)
        vaf_df_sample$subclonal_mutationload <- subclonal_mutationload
      }
      sample_filtered_vcf_object_fraction <- sample_vcf[rownames(sample_vcf) %in% 
                                                           rownames(df_read_based_fraction), ]
      #Make mutation matrix
      fraction_mut_mat <- mut_matrix(list(rowRanges(sample_filtered_vcf_object_fraction)),
                                     ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
      colnames(fraction_mut_mat) <- paste0(anonymous_id, "Dx_", fraction)
      #Add mutation matrix to large matrix
      mut_mat_collection_combined <- cbind(mut_mat_collection_combined, fraction_mut_mat)
      #Save files
      writeVcf(
        sample_filtered_vcf_object_fraction,
        paste0(
          "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/filteredVcf_diagnosis/",
          paste0(anonymous_id, "Dx_", fraction),
          ".vcf"
        )
      )
      write.csv(
        fraction_mut_mat,
        paste0(
          "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/DiagnosisMutationMatrices/",
          paste0(anonymous_id, "Dx_", fraction),
          ".csv"
        )
      )
    }
  } else {
    #Fill in columns
    vaf_df_sample$cutoff <- NA
    vaf_df_sample$subclonal_mutationload <- NA
    vaf_df_sample$clonal_mutationload <- nrow(df_read_based_filtered_all)
    #Make mutation matrix
    sample_filtered_vcf_object <- sample_vcf[rownames(sample_vcf) %in%
                                               rownames(df_read_based_filtered_all), ]
    sample_mut_mat <- mut_matrix(list(rowRanges(sample_filtered_vcf_object)),
                                   ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
    colnames(sample_mut_mat) <- paste0(anonymous_id, "Dx_clonal")
    #Add mutation matrix to large matrix
    mut_mat_collection_combined <- cbind(mut_mat_collection_combined, sample_mut_mat)
  }
  #Add to combined dataframe
  vaf_df <- rbind(vaf_df, vaf_df_sample)
}
vaf_df <- vaf_df[!is.na(vaf_df$Sample),]

#Get mutational matrices
mut_mat_collection_subclonal <- mut_mat_collection_combined[, grep("_subclonal", colnames(mut_mat_collection_combined))]
mut_mat_collection_clonal <- mut_mat_collection_combined[, grep("_clonal", colnames(mut_mat_collection_combined))]
colnames(mut_mat_collection_subclonal) <- paste0(colnames(mut_mat_collection_subclonal),
                                                 "\nn = ",
                                                 colSums(mut_mat_collection_subclonal))
colnames(mut_mat_collection_clonal) <- paste0(colnames(mut_mat_collection_clonal),
                                              "\nn = ",
                                              colSums(mut_mat_collection_clonal))

#Plot vaf
vaf_df$Modality <- factor(vaf_df$Modality, levels = c("bimodal",
                                                      "multimodal",
                                                      "unimodal"))
vaf_df$cutoff_randomlysampled <- unlist(lapply(vaf_df$cutoff, sample, size = 1))
vaf_df <- arrange(vaf_df, Modality, desc(subclonal_mutationload))
vaf_df$SampleFactor <- factor(vaf_df$Sample, levels = unique(vaf_df$Sample))
pdf("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/ViolinPlots.pdf",
    width = 10,
    height = 10)
ggplot(vaf_df, aes(x = "", y = VAF)) +
  geom_violin(fill = "snow3") +
  facet_wrap(~SampleFactor) +
  ylim(0,1) +
  geom_hline(aes(yintercept = cutoff_randomlysampled)) +
  geom_label(mapping = aes(x = -Inf, y = 0.95, label = SampleFactor),
            hjust = -0.1,
            distinct(vaf_df, SampleFactor)) +
  geom_text(mapping = aes(x = Inf, y = 0.95, label = paste0("n=", clonal_mutationload)),
             hjust = 1.1,
             distinct(vaf_df, SampleFactor, clonal_mutationload),
            size = 3.5) +
  geom_text(mapping = aes(x = Inf, y = 0.1, label = ifelse(!is.na(subclonal_mutationload),
                                                            paste0("n=", subclonal_mutationload),
                                                            subclonal_mutationload)),
             hjust = 1.1,
             distinct(vaf_df, SampleFactor, subclonal_mutationload),
            size = 3.5) +
  theme_bw() +
  theme(strip.text.x = element_blank())
dev.off()

#Get signatures
signatures <- read.table("~/Projects/hypermutated_ALL/RESULTS/tables/DeNovoSignatures_SBS7aALLNoClusters.tsv")
signatures_matrix <- as.matrix(signatures)

#Perform refit
fit_res_subclonal <- fit_to_signatures_strict(mut_mat_collection_subclonal,
                                              signatures_matrix,
                                              max_delta = 0.033,
                                              method = "best_subset")
fit_res_subclonal <- fit_res_subclonal$fit_res
plot_contribution(fit_res_subclonal$contribution)
fit_res_clonal <- fit_to_signatures_strict(mut_mat_collection_clonal,
                                           signatures_matrix,
                                           max_delta = 0.033,
                                           method = "best_subset")
fit_res_clonal <- fit_res_clonal$fit_res
plot_contribution(fit_res_clonal$contribution)

SBS7a_contribution <- data.frame(contribution = fit_res_subclonal$contribution["SBS7a",]/colSums(fit_res_subclonal$contribution),
                                 sample = str_remove(names(fit_res_subclonal$contribution["SBS7a",]), "_.+\\n.+"),
                                 vaf = "subclonal",
                                 abs_contribution = fit_res_subclonal$contribution["SBS7a",]/colSums(fit_res_subclonal$contribution)*colSums(mut_mat_collection_subclonal))
SBS7a_contribution <- rbind(SBS7a_contribution,
                            data.frame(contribution = fit_res_clonal$contribution["SBS7a",]/colSums(fit_res_clonal$contribution),
                                       sample = str_remove(names(fit_res_clonal$contribution["SBS7a",]), "_.+\\n.+"),
                                       vaf = "clonal",
                                       abs_contribution = fit_res_clonal$contribution["SBS7a",]/colSums(fit_res_clonal$contribution)*colSums(mut_mat_collection_clonal)))

#Pick a few patients who clearly show subclones
selected_samples <- unique(vaf_df$Sample[vaf_df$Modality == "bimodal" &
                                           vaf_df$subclonal_mutationload >= 200])
vaf_df_selected <- vaf_df[vaf_df$Sample %in% selected_samples,]

#Check mutational profiles
mut_mat_collection_subclonal_selected <- mut_mat_collection_subclonal[, str_remove(colnames(mut_mat_collection_subclonal), "_.+\n.+") %in% selected_samples]
mut_mat_collection_clonal_selected <- mut_mat_collection_clonal[, str_remove(colnames(mut_mat_collection_clonal), "_.+\n.+") %in% selected_samples]

#Perform refit
fit_res_subclonal_selected <- fit_to_signatures_strict(mut_mat_collection_subclonal_selected,
                                                       signatures_matrix,
                                                       max_delta = 0.033,
                                                       method = "best_subset")
fit_res_subclonal_selected <- fit_res_subclonal_selected$fit_res
plot_contribution(fit_res_subclonal_selected$contribution)
fit_res_clonal_selected <- fit_to_signatures_strict(mut_mat_collection_clonal_selected,
                                                    signatures_matrix,
                                                    max_delta = 0.033,
                                                    method = "best_subset")
fit_res_clonal_selected <- fit_res_clonal_selected$fit_res
plot_contribution(fit_res_clonal_selected$contribution)

SBS7a_contribution_selected <- data.frame(contribution = fit_res_subclonal_selected$contribution["SBS7a",]/colSums(fit_res_subclonal_selected$contribution),
                                          sample = str_remove(names(fit_res_subclonal_selected$contribution["SBS7a",]), "_.+\\n.+"),
                                          vaf = "subclonal",
                                          abs_contribution = fit_res_subclonal_selected$contribution["SBS7a",]/colSums(fit_res_subclonal_selected$contribution)*colSums(mut_mat_collection_subclonal_selected))
SBS7a_contribution_selected <- rbind(SBS7a_contribution_selected,
                                     data.frame(contribution = fit_res_clonal_selected$contribution["SBS7a",]/colSums(fit_res_clonal_selected$contribution),
                                                sample = str_remove(names(fit_res_clonal_selected$contribution["SBS7a",]), "_.+\\n.+"),
                                                vaf = "clonal",
                                                abs_contribution = fit_res_clonal_selected$contribution["SBS7a",]/colSums(fit_res_clonal_selected$contribution)*colSums(mut_mat_collection_clonal_selected)))

#Make sample specific plots
make_vaf_overviewplot <- function(sample_id){
  VAF_plot <- ggplot(vaf_df[vaf_df$Sample == sample_id,], aes(x = "", y = VAF)) +
    geom_violin(fill = "snow3") +
    ylim(0,1) +
    geom_hline(aes(yintercept = cutoff_randomlysampled)) +
    geom_label(mapping = aes(x = -Inf, y = 0.95, label = SampleFactor),
               hjust = -0.1,
               distinct(vaf_df[vaf_df$Sample == sample_id,], SampleFactor)) +
    geom_text(mapping = aes(x = Inf, y = 0.95, label = paste0("n=", clonal_mutationload)),
              hjust = 1.1,
              distinct(vaf_df[vaf_df$Sample == sample_id,], SampleFactor, clonal_mutationload),
              size = 3.5) +
    geom_text(mapping = aes(x = Inf, y = 0.1, label = ifelse(!is.na(subclonal_mutationload),
                                                             paste0("n=", subclonal_mutationload),
                                                             subclonal_mutationload)),
              hjust = 1.1,
              distinct(vaf_df[vaf_df$Sample == sample_id,], SampleFactor, subclonal_mutationload),
              size = 3.5) +
    theme_bw() +
    theme(strip.text.x = element_blank(),
          axis.title.y = element_text(size = 8)) +
    xlab("")
  mut_mat_patient <- mut_mat_collection_combined[,grep(sample_id,
                                                       colnames(mut_mat_collection_combined))]
  colnames(mut_mat_patient) <- colnames(mut_mat_patient) %>%
    str_remove(paste0(sample_id, "_")) %>%
    str_remove("_.+")
  profile_plot <- plot_96_profile(mut_mat_patient,
                                  ymax = 0.3, condensed = T) +
    ylab("")
  contribution_plot <- ggplot(SBS7a_contribution[grep(sample_id, SBS7a_contribution$sample),]) +
    geom_col(aes(x = vaf, y = abs_contribution), fill = cbp2[4]) +
    theme_bw() +
    xlab("") +
    ylab("Absolute Contribution SBS7a") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 8)) +
    ylim(0,3000)
  fullplot <- plot_grid(VAF_plot, profile_plot, contribution_plot,
                        ncol = 3, rel_widths = c(2,3,1), align = "tb")
  return(fullplot)
}

fullplot_P0621 <- make_vaf_overviewplot("P0621Dx")
fullplot_P0628 <- make_vaf_overviewplot("P0628Dx")
fullplot_P0633 <- make_vaf_overviewplot("P0633Dx")
fullplot_P0629 <- make_vaf_overviewplot("P0629Dx")
fullplot_P0429 <- make_vaf_overviewplot("P0429Dx")
fullplot_P0557 <- make_vaf_overviewplot("P0557Dx")

#Save plots so they can be combined in a figure
save(object = fullplot_P0621, file = "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0621.rdata")
save(object = fullplot_P0633, file = "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0633.rdata")
save(object = fullplot_P0629, file = "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0629.rdata")
save(object = fullplot_P0429, file = "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0429.rdata")
save(SBS7a_contribution_selected, file = "~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/DiagnosisMutationMatrices/SBS7aContribution_selected.rdata")