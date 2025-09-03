#Load packages
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(tidyverse)
library(readxl)
library(Polychrome)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
P60 = createPalette(60,  cbp2)
signature_palette <- c("#DADADA",
                       "#9BC0CD",
                       "#CE2627",
                       "#F6EB13",
                       "#B12325",
                       "#AA6AAC",
                       "#5C6BC0")
metadata <- read_excel("~/Projects/hypermutated_ALL/DOCS/SBS7aCohort_metadata.xlsx")
known_patients <- metadata$`Anonymous ID`[!is.na(metadata$`Anonymous ID`)]

#Load own samples mutational matrix
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/SBS7aNoClusters_MutationMatrix.rdata")

#Change mutations matrix to mut_mat
mut_mat <- sbs_mutation_matrix_combined

#Load own signatures after NMF at HPC
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/nmf_res_SBS7aALLNoClusters_5signatures.rdata")

#Change signatures to nmf_res (any variable works, but this way you don't need to change anything)
nmf_res = nmf_res

#Rename signatures
colnames(nmf_res$signatures) <- paste0("Signature ", LETTERS[1:ncol(nmf_res$signatures)])
rownames(nmf_res$contribution) <- paste0("Signature ", LETTERS[1:ncol(nmf_res$signatures)])

#Check similarity to COSMIC and rename accordingly
signatures = get_known_signatures()
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.95)
colnames(nmf_res$signatures)

#Refit de novo signatures to COSMIC if not similar to one COSMIC signature
strict_refit <- fit_to_signatures_strict(nmf_res$signatures, signatures,
                                         max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res

#Make matrix with COSMIC signatures and new signatures (if any left)
#De novo signatures that can be explained by 3 or less COSMIC signatures while having reconstruction > 0.9 are replaced by COSMIC signatures
signatures_denovo = cbind(as.matrix(nmf_res$signatures[,"SBSA"]),
                          signatures[,c("SBS1", "SBS2", "SBS7a", "SBS13", "SBS18", "SBS87")])
colnames(signatures_denovo)[1] <- "SBSA"

#Scale signatures so they all add up to 1
signatures_denovo_scaled <- signatures_denovo
for (col in 1:ncol(signatures_denovo)){
  signatures_denovo_scaled[, col] <- signatures_denovo[, col] / sum(signatures_denovo[, col])
}
apply(signatures_denovo_scaled, 2, sum)
write.table(signatures_denovo_scaled, "~/Projects/hypermutated_ALL/RESULTS/tables/DeNovoSignatures_SBS7aALLNoClusters.tsv", sep = "\t")

#Refit own samples with matrix COSMIC/new signatures
strict_refit <- fit_to_signatures_strict(mut_mat, signatures_denovo_scaled, max_delta = 0.033, method = "best_subset")
fit_res_strict <- strict_refit$fit_res
fit_res_strict$signatures <- signatures_denovo_scaled
save(object = fit_res_strict, file = "~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/DeNovo_Refit_SBS7aALLNoClusters.rdata")

#Get samples that have a reconstruction value of at least 0.9
good_reconstructions <- rownames(as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))))[as.data.frame(round(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed)), digits = 2)) >= 0.9]

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
write.csv(relative_contributions, "~/Projects/hypermutated_ALL/RESULTS/tables/RelativeContributions_SBS7aALLNoClusters.csv")
write.csv(relative_contributions[, good_reconstructions], "~/Projects/hypermutated_ALL/RESULTS/tables/RelativeContributions_SBS7aALLNoClusters_goodreconstruction.csv")
SBS7_samples_notstrict <- colnames(relative_contributions[,relative_contributions["SBS7a",]>0])
SBS7_samples_notstrict <- SBS7_samples_notstrict[order(SBS7_samples_notstrict)]

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
write.csv(absolute_contributions, "~/Projects/hypermutated_ALL/RESULTS/tables/AbsoluteContributions_SBS7aALLNoClusters.csv")
absolute_contributions[,!grepl("P_", colnames(absolute_contributions))]

#Perform bootstrapping
contri_boots <- fit_to_signatures_bootstrapped(mut_mat[, SBS7_samples_notstrict],
                                               signatures_denovo_scaled,
                                               n_boots = 100,
                                               method = "strict_best_subset",
                                               max_delta = 0.033)
plot_bootstrapped_contribution(contri_boots) +
  facet_grid(sample ~ ., scales = "free")

#Make final list of SBS7a samples, which have a reconstruction value of at least 0.9 and bootstrap of at least 75%
SBS7_samples_reconstruction <- SBS7_samples_notstrict[SBS7_samples_notstrict %in% good_reconstructions]
bootstrap_overview <- table(str_remove(rownames(contri_boots), "_\\d+$")[contri_boots[, "SBS7a"] == 0])
SBS7_samples_complete <- SBS7_samples_reconstruction[!SBS7_samples_reconstruction %in% names(bootstrap_overview)[bootstrap_overview > 10]]