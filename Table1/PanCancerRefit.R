#Load packages
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(tidyverse)
library(readxl)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
library(Polychrome)
P60 = createPalette(60,  cbp2)
P12 = c("grey", createPalette(11,  cbp2))

#Load mutational matrix of all samples
load("~/Projects/hypermutated_ALL/ANALYSES/biobankData/PanCancer/Data_Maxima_WGS/pancancerbiobank_mut_mat_05082024.rdata")

#Load mutational matrix of filtered samples
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/SBS7a_pancancerbiobank_mut_mat28042025.rdata")

#Load metadata
metadata <- read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/Lianne/UV/UVpatientsOverview.xlsx")
UV_bms <- metadata$`BioID Tumor`
UV_bms <- UV_bms[!is.na(UV_bms)]
UV_bms <- unlist(str_split(UV_bms, ","))

#Get biobank version of UV samples
sbs_mutation_matrix_combined <- sbs_mutation_matrix_combined[,
                                                             grep("^PM",
                                                                  colnames(sbs_mutation_matrix_combined))]
mut_mat <- cbind(sbs_mutation_matrix_combined,
                 pancancer_mutation_matrix[, str_remove(colnames(pancancer_mutation_matrix), "_.+") %in%
                                             UV_bms])

#Load own signatures after NMF at HPC
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/nmf_res_SBS7a_pancancerbiobank_28042025_8signatures.rdata")

#Change signatures to nmf_res (any variable works, but this way you don't need to change anything)
nmf_res = nmf_res

#Rename signatures to some name
colnames(nmf_res$signatures) <- paste0("Signature ", LETTERS[1:ncol(nmf_res$signatures)])
rownames(nmf_res$contribution) <- paste0("Signature ", LETTERS[1:ncol(nmf_res$signatures)])

#Check similarity to COSMIC and rename accordingly
signatures = get_known_signatures()
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.95)
colnames(nmf_res$signatures)

#Plot signatures
plot_96_profile(nmf_res$signatures, ymax = 0.3)

#Refit de novo signatures to COSMIC if not similar to one COSMIC signature
strict_refit <- fit_to_signatures_strict(nmf_res$signatures, signatures, max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res

#See signatures in refit plotted in barplots
plot_contribution(fit_res_strict$contribution[apply(fit_res_strict$contribution, 1, sum) != 0,],
                  coord_flip = FALSE,
                  mode = "relative",
                  palette = as.vector(P60))

#Check cosine similarity between de novo signatures and reconstructed profile from refit with COSMIC
plot_original_vs_reconstructed(nmf_res$signatures[], fit_res_strict$reconstructed[], 
                               y_intercept = 0.9)

#Make matrix with COSMIC signatures and new signatures (if any left)
signatures_denovo = cbind(as.matrix(nmf_res$signatures[,"SBSE"]),
                          signatures[,c("SBS1", "SBS2", "SBS5", "SBS7a", "SBS8", "SBS9",
                                        "SBS13", "SBS15", "SBS17b", "SBS18")])
colnames(signatures_denovo)[1] <- "SBSE"
signatures_denovo_scaled <- signatures_denovo
for (col in 1:ncol(signatures_denovo)){
  signatures_denovo_scaled[, col] <- signatures_denovo[, col] / sum(signatures_denovo[, col])
}
apply(signatures_denovo_scaled, 2, sum)

#Refit own samples with matrix COSMIC/new signatures
strict_refit <- fit_to_signatures_strict(mut_mat, signatures_denovo_scaled, max_delta = 0.033, method = "best_subset")
fit_res_strict <- strict_refit$fit_res
fit_res_strict$signatures <- signatures_denovo_scaled
save(object = fit_res_strict, file = "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/DeNovo_Refit_SBS7a_pancancerbiobank_8signatures_28042025.rdata")

#Calculate relative contribution
contribution_with_totals <- rbind(fit_res_strict$contribution, apply(fit_res_strict$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col)-1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

#Calculate absolute contributions
contribution_with_original_totals <- rbind(contribution_with_totals, apply(mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"
calculate_absolute_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col)-2)]
  refit_total <- contribution_col[length(contribution_col)-1]
  real_total <- contribution_col[length(contribution_col)]
  absolute_contributions <- (contributions / refit_total) * real_total
  return(absolute_contributions)
}
absolute_contributions <- apply(contribution_with_original_totals, 2, calculate_absolute_contribution)

#Select samples with some SBS7a contribution
SBS7_samples_notstrict <- colnames(relative_contributions[,relative_contributions["SBS7a",]>0])
SBS7_samples_notstrict <- SBS7_samples_notstrict[order(SBS7_samples_notstrict)]

#Perform bootstrapping
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[, SBS7_samples_notstrict],
                                               signatures_denovo_scaled,
                                               n_boots = 100,
                                               method = "strict_best_subset",
                                               max_delta = 0.033)
plot_bootstrapped_contribution(contri_boots) +
  facet_grid(sample ~ ., scales = "free")

#Save bootstrap results
save(object = contri_boots,
     file = "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/Bootstrap_SBS7a_pancancerbiobank_8signatures_28042025.rdata")

#Get final list of samples
good_reconstructions <- rownames(as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))))[as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))) >= 0.9]
SBS7_samples_reconstruction <- SBS7_samples_notstrict[SBS7_samples_notstrict %in% good_reconstructions]
bootstrap_overview <- table(str_remove(rownames(contri_boots), "_\\d+$")[contri_boots[, "SBS7a"] == 0])
plot_bootstrapped_contribution(contri_boots[str_remove(rownames(contri_boots), "_\\d+$") %in% 
                                              names(bootstrap_overview)[bootstrap_overview > 10], ]) +
  facet_grid(sample ~ ., scales = "free")
SBS7_samples_complete <- SBS7_samples_reconstruction[!SBS7_samples_reconstruction %in% names(bootstrap_overview)[bootstrap_overview > 10]]
str_extract(SBS7_samples_complete, "\\d+") %>% unique() %>% length()

#Calculate percentage of SBS7a-positive patients per tumor type
tumor_type_counts <- as.data.frame(table(str_remove((str_remove(colnames(pancancer_mutation_matrix), ".+\\(")), "\\)")))
colnames(tumor_type_counts) <- c("Tumor Type", "Number of Samples")
tumor_type_counts <- tumor_type_counts[!tumor_type_counts$`Tumor Type` %in% c("Benign", "MDS", "Other", "Uncertain", "B-cell unspecified", "Leukemia unspecified", "Lymphoma unspecified"),]
tumor_type_counts$SBS7a_samples <- 0
for (tumor_type in tumor_type_counts$`Tumor Type`){
  if (tumor_type == "BALL"){
    ALL_count1 <- sum(unique(str_extract(grep("-", SBS7_samples_complete, value = T, invert = T), "\\d+")) %in% str_extract(colnames(pancancer_mutation_matrix), "\\d\\d\\d\\d+"))
    ALL_count2 <- table(str_remove(str_remove(SBS7_samples_complete, ".+\\("), "\\)"))[tumor_type]
    tumor_type_counts$SBS7a_samples[tumor_type_counts$`Tumor Type` == tumor_type] <- ALL_count1 + ALL_count2
  } else if (tumor_type %in% str_remove(str_remove(SBS7_samples_complete, ".+\\("), "\\)")){
    tumor_type_counts$SBS7a_samples[tumor_type_counts$`Tumor Type` == tumor_type] <- table(str_remove(str_remove(SBS7_samples_complete, ".+\\("), "\\)"))[tumor_type]
  }
}
tumor_type_counts$SBS7a_percentage <- tumor_type_counts$SBS7a_samples / tumor_type_counts$`Number of Samples` * 100

#Get percentage of SBS7a-negative patients
tumor_type_counts$SBS7anegative_samples <- tumor_type_counts$`Number of Samples` - tumor_type_counts$SBS7a_samples
tumor_type_counts_long <- pivot_longer(tumor_type_counts, c("SBS7a_samples", "SBS7anegative_samples"), names_to = "SBS7a Status", values_to = "Number of patients")
tumor_type_counts_long$`SBS7a Status`[tumor_type_counts_long$`SBS7a Status` == "SBS7a_samples"] <- "SBS7a-positive"
tumor_type_counts_long$`SBS7a Status`[tumor_type_counts_long$`SBS7a Status` == "SBS7anegative_samples"] <- "SBS7a-negative"
tumor_type_counts_long <- arrange(tumor_type_counts_long, desc(`Number of Samples`))
tumor_type_counts_long$`Tumor Type` <- factor(tumor_type_counts_long$`Tumor Type`, levels = unique(tumor_type_counts_long$`Tumor Type`))

#Save source data
write_tsv(tumor_type_counts, "~/Projects/hypermutated_ALL/ANALYSES/biobankData/PanCancer/SBS7apositivePercentage.tsv")

#Classify patients into broader tumor groups
hematological <- c("AML", "B-ALL", "B-cell lymphoma", "CML", "Hodgkin lymphoma",
                   "T-ALL", "T-cell lymphoma")
neurological <- c("Astrocytoma", "Astrocytoma high-grade", "ATRT", "Brain high-grade",
                  "Brain other", "Ependymoma", "Germ cell brain", "Glioma",
                  "Glioma high-grade", "Medulloblastoma")
solid <- c("Adrenocortical carcinoma", "Carcinoma other", "Ewing sarcoma",
           "Fibrosarcoma", "Germ cell other", "Gonadal carcinoma", "Gonadal germ cell",
           "Gonadal other", "Hepatic carcinoma", "Hepatoblastoma", "Histiocytic neoplasm",
           "Melanoma", "Nasopharyngeal carcinoma", "Nephroblastoma", "Neuroblastoma",
           "Osteosarcoma", "Renal carcinoma", "Rhabdomyosarcoma", "Sarcoma other",
           "Thyroid carcinoma")

tumor_type_counts$`Tumor Group`[tumor_type_counts$`Tumor Type` %in% hematological] <- "Hematological Tumors"
tumor_type_counts$`Tumor Group`[tumor_type_counts$`Tumor Type` %in% neurological] <- "Neurological Tumors"
tumor_type_counts$`Tumor Group`[tumor_type_counts$`Tumor Type` %in% solid] <- "Solid Tumors"

#Elongate table
tumor_type_counts_long <- pivot_longer(tumor_type_counts, c("SBS7a_samples", "SBS7anegative_samples"), names_to = "SBS7a Status", values_to = "Number of patients")
tumor_type_counts_long$`SBS7a Status`[tumor_type_counts_long$`SBS7a Status` == "SBS7a_samples"] <- "SBS7a-positive"
tumor_type_counts_long$`SBS7a Status`[tumor_type_counts_long$`SBS7a Status` == "SBS7anegative_samples"] <- "SBS7a-negative"
tumor_type_counts_long <- arrange(tumor_type_counts_long, desc(`Number of Samples`))
tumor_type_counts_long$`Tumor Type` <- factor(tumor_type_counts_long$`Tumor Type`, levels = unique(tumor_type_counts_long$`Tumor Type`))

#Make dataframe per tumor group
tumor_group_counts <- data.frame("Tumor Group" = c("Hematological Tumors",
                                                   "Neurological Tumors",
                                                   "Solid Tumors"),
                                 "SBS7a_samples" = c(sum(tumor_type_counts$SBS7a_samples[tumor_type_counts$`Tumor Group` == "Hematological Tumors"]),
                                                     sum(tumor_type_counts$SBS7a_samples[tumor_type_counts$`Tumor Group` == "Neurological Tumors"]),
                                                     sum(tumor_type_counts$SBS7a_samples[tumor_type_counts$`Tumor Group` == "Solid Tumors"])),
                                 "SBS7anegative_samples" = c(sum(tumor_type_counts$SBS7anegative_samples[tumor_type_counts$`Tumor Group` == "Hematological Tumors"]),
                                                             sum(tumor_type_counts$SBS7anegative_samples[tumor_type_counts$`Tumor Group` == "Neurological Tumors"]),
                                                             sum(tumor_type_counts$SBS7anegative_samples[tumor_type_counts$`Tumor Group` == "Solid Tumors"])),
                                 "Number of Samples" = c(sum(tumor_type_counts$`Number of Samples`[tumor_type_counts$`Tumor Group` == "Hematological Tumors"]),
                                                         sum(tumor_type_counts$`Number of Samples`[tumor_type_counts$`Tumor Group` == "Neurological Tumors"]),
                                                         sum(tumor_type_counts$`Number of Samples`[tumor_type_counts$`Tumor Group` == "Solid Tumors"])))
colnames(tumor_group_counts) <- c("Tumor Group",
                                  "SBS7a_samples",
                                  "SBS7anegative_samples",
                                  "Number of Samples")
tumor_group_counts$SBS7a_percentage <- tumor_group_counts$SBS7a_samples / tumor_group_counts$`Number of Samples` * 100

tumor_group_counts_long <- pivot_longer(tumor_group_counts, c("SBS7a_samples", "SBS7anegative_samples"), names_to = "SBS7a Status", values_to = "Number of patients")
tumor_group_counts_long$`SBS7a Status`[tumor_group_counts_long$`SBS7a Status` == "SBS7a_samples"] <- "SBS7a-positive"
tumor_group_counts_long$`SBS7a Status`[tumor_group_counts_long$`SBS7a Status` == "SBS7anegative_samples"] <- "SBS7a-negative"
tumor_group_counts_long <- arrange(tumor_group_counts_long, desc(`Number of Samples`))
tumor_group_counts_long$`Tumor Group` <- factor(tumor_group_counts_long$`Tumor Group`, levels = unique(tumor_group_counts_long$`Tumor Group`))

#Change table a bit
for (tumor_group in unique(tumor_type_counts$`Tumor Group`)){
  number_of_group_samples <- sum(tumor_type_counts$`Number of Samples`[tumor_type_counts$`Tumor Group` == tumor_group])
  if (tumor_group == unique(tumor_type_counts$`Tumor Group`)[1]){
    tumor_group_counts_short <- data.frame("Tumor Group" = tumor_group,
                                           "Number of Samples" = number_of_group_samples)
  } else {
    tumor_group_counts_short <- rbind(tumor_group_counts_short,
                                      c(tumor_group, number_of_group_samples))
  }
}
colnames(tumor_group_counts_short) <- c("Tumor Group", "Number of Samples")
write_tsv(tumor_group_counts_short, "~/Projects/adultCancerPredisposingGenes/DOCS/NegativeControls.tsv")
