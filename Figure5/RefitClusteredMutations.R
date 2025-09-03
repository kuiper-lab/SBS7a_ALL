#Load packages
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(tidyverse)
library(readxl)
library(Polychrome)

## load own samples mutational matrix
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/pmc_mut_mat_clusters_nogermlinefilter_nosoftclippedbases_015_23042025.rdata")

##change mutations matrix to mut_mat
mut_mat <- sbs_mutation_matrix_pmc

#Load signatures
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/DeNovo_Refit_SBS7aALLNoClusters.rdata")
signatures_denovo_scaled <- fit_res_strict$signatures

## refit own samples with matrix COSMIC/new signatures
strict_refit <- fit_to_signatures_strict(mut_mat, signatures_denovo_scaled, max_delta = 0.033, method = "best_subset")
fit_res_strict <- strict_refit$fit_res
fit_res_strict$signatures <- signatures_denovo_scaled
save(object = fit_res_strict, file = "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/DeNovo_Refit_SBS7aALLClusters_ReuseSignatures.rdata")

## check cosine similarity between original profiles and reconstruction from refit
good_reconstructions <- rownames(as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))))[as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))) >= 0.9]

## calculate relative contribution
contribution_with_totals <- rbind(fit_res_strict$contribution, apply(fit_res_strict$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col)-1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)
write.csv(relative_contributions, "~/Projects/hypermutated_ALL/RESULTS/tables/RelativeContributions_SBS7aALLClusters_ReuseSignatures.csv")
write.csv(relative_contributions[, good_reconstructions], "~/Projects/hypermutated_ALL/RESULTS/tables/RelativeContributions_SBS7aALLClusters_goodreconstruction_ReuseSignatures.csv")
SBS7_samples_notstrict <- colnames(relative_contributions[,relative_contributions["SBS7a",]>0])
SBS7_samples_notstrict <- SBS7_samples_notstrict[order(SBS7_samples_notstrict)]

#Calculate absolute contribution
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
write.csv(absolute_contributions, "~/Projects/hypermutated_ALL/RESULTS/tables/AbsoluteContributions_SBS7aALLClusters_ReuseSignatures.csv")

#Perform bootstrapping
SBS7_samples_reconstruction <- SBS7_samples_notstrict[SBS7_samples_notstrict %in% good_reconstructions]
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[, SBS7_samples_reconstruction],
                                               signatures_denovo_scaled,
                                               n_boots = 100,
                                               method = "strict_best_subset",
                                               max_delta = 0.033)
plot_bootstrapped_contribution(contri_boots) +
  facet_grid(sample ~ ., scales = "free")

#Make final list of SBS7a samples
SBS7_samples_reconstruction <- SBS7_samples_notstrict[SBS7_samples_notstrict %in% good_reconstructions]
bootstrap_overview <- table(str_remove(rownames(contri_boots), "_\\d+$")[contri_boots[, "SBS7a"] == 0])
SBS7_samples_complete <- SBS7_samples_reconstruction[!SBS7_samples_reconstruction %in% names(bootstrap_overview)[bootstrap_overview > 10]]
str_extract(SBS7_samples_complete, "\\d+") %>% unique() %>% length()