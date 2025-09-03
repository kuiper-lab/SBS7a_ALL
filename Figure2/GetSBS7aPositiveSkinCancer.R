#Load packages
library(tidyverse)
library(MutationalPatterns)
library(ggsignif)
library(readxl)

#Load data
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_mut_mat.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_fit_res_strict.rdata")

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

#Select samples with some SBS7a contribution
SBS7_samples_notstrict <- colnames(relative_contributions[,relative_contributions["SBS7a",]>0])
save(object = SBS7_samples_notstrict,
     file = "~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_SBS7_samples_notstrict.rdata")

#Load bootstraps (these were run on the cluster with the same settings as for the BCP-ALL cohort)
bootstrap_files <- list.files("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/bootstrap_ALL_ALCL_skincancer_29042025/",
                              full.names = T)
signatures <- rownames(relative_contributions)
for (bootstrap_file in bootstrap_files){
  load(bootstrap_file)
  #Add absent signatures if necessary
  if (sum(!signatures %in% colnames(contri_boots)) > 0){
    missing_signatures <- signatures[!signatures %in% colnames(contri_boots)]
    empty_signatures <- matrix(rep(0, nrow(contri_boots)*length(missing_signatures)),
                                   ncol = length(missing_signatures), dimnames = list(rownames(contri_boots), missing_signatures))
    contri_boots <- cbind(contri_boots,
                          empty_signatures)
  }
  contri_boots <- contri_boots[, sort(colnames(contri_boots))]
  if (bootstrap_file == bootstrap_files[1]){
    contri_boots_combined <- contri_boots
  } else {
    contri_boots_combined <- rbind(contri_boots_combined,
                                   contri_boots)
  }
}
contri_boots <- contri_boots_combined

#Get final list of samples
good_reconstructions <- rownames(as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))))[as.data.frame(diag(cos_sim_matrix(mut_mat, fit_res_strict$reconstructed))) >= 0.9]
SBS7_samples_reconstruction <- SBS7_samples_notstrict[SBS7_samples_notstrict %in% good_reconstructions]
bootstrap_overview <- table(str_remove(rownames(contri_boots), "_\\d+$")[contri_boots[, "SBS7a"] == 0])
SBS7_samples_complete <- SBS7_samples_reconstruction[!SBS7_samples_reconstruction %in% names(bootstrap_overview)[bootstrap_overview > 10]]
str_extract(SBS7_samples_complete, "\\d+") %>% unique() %>% length()
save(object = SBS7_samples_complete,
     file = "~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_SBS7_samples.rdata")

