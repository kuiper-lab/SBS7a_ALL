#Load packages -----------------------------------------------------------------
library(tidyverse)
library(ggtree)
library(treeio)
library(stringi)
library(stringdist)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggnewscale)
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38'
library(MutationalPatterns)
library(VariantAnnotation)
theme_set(theme_classic())
library(cellPhyWrapperPlotting)
library(bedr)
library(cowplot)
library(ggsignif)
library(car)
source("fit_to_signatures_strict_tree_method.R") #Present in same directory as this script
source("plot_tree_contributions_branchnames.R") #Present in same directory as this script

#Make tree plot ----------------------------------------------------------------
#Load tree and vcf file
tree = readRDS("~/Projects/SingleCellWGS/DATA/CellPhyWrapper/P0625_includingbulktumor_MQandadjustedVAF025filtered_blacklistfiltered/TreeObject0.4.RDS")
vcf = VariantAnnotation::readVcf("~/Projects/SingleCellWGS/DATA/includebulktumor/P0625/MQandadjustedVAF025filtered/GenoType/P0625.output.snv.blacklistfiltered.vcf")

#Prepare tree
tree = prepare_tree(tree)

#Get VCF per branch
branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)

#Get profile per branch
branch_mm = mut_matrix(branch_grl, ref_genome)

#Get contributions per branch
signatures <- read.table("~/Projects/hypermutated_ALL/RESULTS/tables/DeNovoSignatures_SBS7aALLNoClusters.tsv")
signatures_matrix <- as.matrix(signatures)
pta_sig <- read.table("~/Projects/SingleCellWGS/REFERENCES/PTA_signature/PTA_Artefact_Signature.txt",
                      sep = "\t", header = T, row.names = 1)
pta_sig = as.matrix(pta_sig)
sigs_slct = cbind(signatures_matrix, pta_sig)
contribution = fit_to_signatures_strict_tree(mut_matrix = branch_mm,
                                             signatures = sigs_slct,
                                             max_delta = 0.033,
                                             remove_min = 100,
                                             method = "best_subset")

#Add signature contributions to your tree
tree = add_contribution(tree, contribution = contribution)

# rename the branches
old = tree@phylo$tip.label # get old labels
new = str_replace(old, "P0625", "P0625")
new[!grepl("P0625", new)] <- "Tumor Bulk"
tip_names_replace = data.frame(old, new)
tree_rename = change_tip_labels(tree = tree, id_change_tb = tip_names_replace)

#Plot
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625_Tree_SBS7a_MQandadjustedVAF025filtered_blacklistfiltered_includingbranchnames.pdf",
    width = 10)
plot_tree_contribution_branchnames(tree_rename, common_name = 'Lineage Tree',
                                   add_branch_length = T,
                                   add_bootstrap = F, signature = 'SBS7a', pie_size = 1,
                                   branch_names = T,
                                   remove_min = 100)
dev.off()
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625_Tree_SBS1_MQandadjustedVAF025filtered_blacklistfiltered_includingbranchnames.pdf",
    width = 10)
plot_tree_contribution_branchnames(tree_rename, common_name = 'Lineage Tree',
                                   add_branch_length = T,
                                   add_bootstrap = F, signature = 'SBS1', pie_size = 1,
                                   branch_names = T,
                                   remove_min = 100)
dev.off()
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625_Tree_PTA_MQandadjustedVAF025filtered_blacklistfiltered_includingbranchnames.pdf",
    width = 10)
plot_tree_contribution_branchnames(tree_rename, common_name = 'Lineage Tree',
                                   add_branch_length = T,
                                   add_bootstrap = F, signature = 'PTA', pie_size = 1,
                                   branch_names = T,
                                   remove_min = 100)
dev.off()

#Plot quality measures ---------------------------------------------------------
##Analyse depth
#Check DP
for (idx in 1:length(branch_vcf)){
  branch <- names(branch_vcf)[idx]
  branchspecific_vcf <- branch_vcf[[idx]]
  #Get branch-specific DP
  branch_DP <- as.data.frame(geno(branchspecific_vcf)$DP)
  if (nrow(branch_DP) > 0){
    #Add some information
    branch_DP$mutation <- rownames(branch_DP)
    branch_DP$branch <- branch
    #Make long
    branch_DP_long <- pivot_longer(branch_DP,
                                   cols = 1:ncol(branchspecific_vcf),
                                   names_to = "Cell",
                                   values_to = "DP")
    #Make combined dataframe
    if (idx == 1){
      DP_df <- branch_DP_long
    } else {
      DP_df <- rbind(DP_df, branch_DP_long)
    }
  }
}

#Get branches that have at least 100 mutations
DP_df %>%
  group_by(branch) %>%
  summarise(count = length(unique(mutation))) ->
  branch_sizes
DP_df_large_branches <- DP_df[DP_df$branch %in%
                                branch_sizes$branch[branch_sizes$count >= 100],]

#Remove bulk samples because they have different depths
DP_df_large_branches_nobulk <- DP_df_large_branches[DP_df_large_branches$Cell != "P0625MSC" & DP_df_large_branches$Cell != "TumorBulk",]

#Put branches in order from left to right
DP_df_large_branches_nobulk$branch <- factor(DP_df_large_branches_nobulk$branch,
                                             levels = c("X",
                                                        "D",
                                                        "F",
                                                        "C",
                                                        "G",
                                                        "R",
                                                        "V",
                                                        "Q",
                                                        "B",
                                                        "A",
                                                        "P",
                                                        "H",
                                                        "J",
                                                        "I",
                                                        "K"))

#Get some statistics
group_by(DP_df_large_branches_nobulk,
         branch) %>%
  summarise(avg = mean(DP)) %>%
  ggplot(aes(x = branch, y = avg)) +
  geom_col() +
  ylab("Mean Depth Across All Cells") +
  xlab("Branch") +
  theme_bw() ->
  mean_depth_plot
group_by(DP_df_large_branches_nobulk,
         branch) %>%
  summarise(median = median(DP)) %>%
  ggplot(aes(x = branch, y = median)) +
  geom_col() +
  ylab("Median Depth Across All Cells") +
  xlab("Branch") +
  theme_bw() ->
  median_depth_plot
group_by(DP_df_large_branches_nobulk,
         branch, Cell) %>%
  summarise(avg = median(DP)) %>%
  group_by(branch) %>%
  summarise(sd = sd(avg)) %>%
  ggplot(aes(x = branch, y = sd)) +
  geom_col() +
  ylab("SD Of Median Depth Per Cell") +
  xlab("Branch") +
  theme_bw() ->
  sd_mediandepthpercell_plot
group_by(DP_df_large_branches_nobulk,
         branch, Cell) %>%
  summarise(avg = mean(DP)) %>%
  group_by(branch) %>%
  summarise(sd = sd(avg)) %>%
  ggplot(aes(x = branch, y = sd)) +
  geom_col() +
  ylab("SD Of Mean Depth Per Cell") +
  xlab("Branch") +
  theme_bw() ->
  sd_meandepthpercell_plot
group_by(DP_df_large_branches_nobulk,
         branch, Cell) %>%
  summarise(avg = median(DP)) %>%
  group_by(branch) %>%
  summarise(min = min(avg)) %>%
  ggplot(aes(x = branch, y = min)) +
  geom_col() +
  ylab("Minimal Median Depth Per Cell") +
  xlab("Branch") +
  theme_bw() ->
  min_mediandepthpercell_plot
group_by(DP_df_large_branches_nobulk,
         branch, Cell) %>%
  summarise(avg = mean(DP)) %>%
  group_by(branch) %>%
  summarise(min = min(avg)) %>%
  ggplot(aes(x = branch, y = min)) +
  geom_col() +
  ylab("Minimal Mean Depth Per Cell") +
  xlab("Branch") +
  theme_bw() ->
  min_meandepthpercell_plot

#Combine into one plot
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625_Quality.pdf",
    width = 10, height = 8)
plot_grid(mean_depth_plot, min_meandepthpercell_plot, sd_meandepthpercell_plot,
          median_depth_plot, min_mediandepthpercell_plot, sd_mediandepthpercell_plot,
          ncol = 2, byrow = F)
dev.off()

#Plot tumor VAF of SBS7a-positive branches -------------------------------------
#Get branches with a contribution of SBS7a
SBS7aPositiveBranches <- colnames(contribution)[contribution["SBS7a",]>0]
for (branch_idx in 1:length(SBS7aPositiveBranches)){
  branch <- SBS7aPositiveBranches[branch_idx]
  branch_specific_vcf <- branch_vcf[[branch]]
  bulk_AD <- geno(branch_specific_vcf)$AD[, "TumorBulk"]
  bulk_ref <- as.numeric(lapply(bulk_AD, "[[", 1))
  bulk_alt <- as.numeric(lapply(bulk_AD, "[[", 2))
  bulk_vaf_branch <- data.frame(VAF = bulk_alt/(bulk_ref + bulk_alt),
                                Branch = branch,
                                chromosome = str_remove(rownames(branch_specific_vcf), ":.+"))
  if (branch_idx == 1){
    bulk_vaf_df <- bulk_vaf_branch
  } else {
    bulk_vaf_df <- rbind(bulk_vaf_df,
                         bulk_vaf_branch)
  }
}

bulk_vaf_df$Branch <- factor(bulk_vaf_df$Branch,
                             levels = c("X", "F", "G", "P"))

pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625SBS7aPositiveBranchesBulk.pdf",
    height = 4)
ggplot(bulk_vaf_df) +
  geom_density(aes(x = VAF, col = Branch)) +
  xlim(0,1) +
  xlab("Allele Frequency In Tumor")
dev.off()

#Add CNA information
bulk_vaf_df$CNA <- NA
bulk_vaf_df$CNA[bulk_vaf_df$chromosome == "1"] <- "p-gain" 
bulk_vaf_df$CNA[bulk_vaf_df$chromosome == "6"] <- "subclonal_gain" 
bulk_vaf_df$CNA[bulk_vaf_df$chromosome == "9"] <- "small_gain" 
bulk_vaf_df$CNA[bulk_vaf_df$chromosome %in% c("2", "13", "16")] <- "LOH"
bulk_vaf_df$CNA[bulk_vaf_df$chromosome %in% c("4", "7", "8", "10", "14", "18", "21")] <- "gain"
bulk_vaf_df$CNA[bulk_vaf_df$chromosome %in% c("3", "5", "11", "12", "15", "17", "19", "20", "22")] <- "noCNA"
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625SBS7aPositiveBranchesBulk_CNAannotated.pdf",
    height = 4)
ggplot(bulk_vaf_df[!is.na(bulk_vaf_df$CNA),]) +
  geom_density(aes(x = VAF, col = Branch)) +
  facet_wrap(~ CNA) +
  xlim(0,1) +
  xlab("Allele Frequency In Tumor")
dev.off()

#Perform tests for mutations that do not overlap with CNAs
bulk_vaf_df_noCNA <- bulk_vaf_df[bulk_vaf_df$CNA == "noCNA" &
                                   !is.na(bulk_vaf_df$CNA),]
#Test for equality of variances and normality of data
leveneTest(VAF ~ Branch, data = bulk_vaf_df_noCNA)
res_aov <- aov(VAF ~ Branch,
               data = bulk_vaf_df_noCNA)
qqPlot(res_aov$residuals,
       id = FALSE)
shapiro.test(res_aov$residuals)
#Variances are not equal so we need to use nonparametric tests
kruskal.test(VAF ~ Branch, bulk_vaf_df_noCNA)
pairwise.wilcox.test(bulk_vaf_df_noCNA$VAF, bulk_vaf_df_noCNA$Branch,
                     p.adjust.method = "bonferroni")
pairwise.wilcox.test(bulk_vaf_df$VAF, bulk_vaf_df$Branch,
                     p.adjust.method = "bonferroni")
#Plot
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625SBS7aPositiveBranchesBulk_boxplot.pdf",
    height = 4, width = 5)
ggplot(bulk_vaf_df_noCNA,
       aes(x = Branch, y = VAF, fill = Branch)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("F", "X"),
                                 c("G", "X"),
                                 c("P", "X")), map_signif_level = T, y_position = c(0.8, 0.86, 0.92)) +
  ylab("Allele Frequency In Tumor") +
  ylim(0,1)
dev.off()

#Plot mutational profiles for branches
pdf("~/Projects/SingleCellWGS/RESULTS/SBS7a/Figures/P0625SBS7aPositiveBranchesMutationalProfile.pdf",
    height = 4)
plot_96_profile(branch_mm[,c("X",
                             "F",
                             "G",
                             "P")],
                condensed = T)
dev.off()
