#Load packages -----------------------------------------------------------------
library(data.table)
library(scales)
library(tidyverse)
library(cowplot)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(readxl)
library(karyoploteR)
library(ggsignif)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


#Read in reference data --------------------------------------------------------
centro <- data.frame(fread('~/Projects/hypermutated_ALL/CODE/R/centromeres/cytoBand.hg38.centromeresOnly.txt'))
chromLength <- data.frame(fread('~/Downloads/chromosomeSizes.hg38.noCommas.txt'))
anonymous_names <- read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx")
metadata <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aPositiveiAMP21.xlsx")
metadata_biobank <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aNegativeiAMP21.xlsx")

#Read in St Jude data ------------------------------------------------------------------
#SBS7a positive cases
SBS7apositive_iamp21_path_stjude <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/StJude/iamp21_SBS7aPositive/"
SBS7apositive_iamp21s_names_stjude <- c()
sample_idx <- 1
for (sample_path in list.files(SBS7apositive_iamp21_path_stjude, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), ".rdata")
  SBS7apositive_iamp21s_names_stjude <- c(SBS7apositive_iamp21s_names_stjude, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7apositive_iamp21s_stjude_gr <- sample_gr
  } else {
    SBS7apositive_iamp21s_stjude_gr <- c(SBS7apositive_iamp21s_stjude_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#SBS7a negative cases
SBS7anegative_iamp21_path_stjude <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/StJude/iamp21_SBS7aNegative/"
SBS7anegative_iamp21s_names_stjude <- c()

sample_idx <- 1
for (sample_path in list.files(SBS7anegative_iamp21_path_stjude, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), ".rdata")
  SBS7anegative_iamp21s_names_stjude <- c(SBS7anegative_iamp21s_names_stjude, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7anegative_iamp21s_stjude_gr <- sample_gr
  } else {
    SBS7anegative_iamp21s_stjude_gr <- c(SBS7anegative_iamp21s_stjude_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#Read in PMC data --------------------------------------------------------------
#SBS7a positive cases
SBS7apositive_iamp21_path_pmc <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/SBS7aPositiveiAMP21s/"
SBS7apositive_iamp21s_names_pmc <- c()
recode_list <- metadata$Anonymous_ID
names(recode_list) <- metadata$Tumor_ID
sample_idx <- 1
for (sample_path in list.files(SBS7apositive_iamp21_path_pmc, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), "_WGS.rdata")
  sample_name <- recode(sample_name,!!!recode_list)
  SBS7apositive_iamp21s_names_pmc <- c(SBS7apositive_iamp21s_names_pmc, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7apositive_iamp21s_pmc_gr <- sample_gr
  } else {
    SBS7apositive_iamp21s_pmc_gr <- c(SBS7apositive_iamp21s_pmc_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#SBS7a negative cases
SBS7anegative_iamp21_path_pmc <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/SBS7aNegativeiAMP21s/"
SBS7anegative_iamp21s_names_pmc <- c()
sample_idx <- 1
for (sample_path in list.files(SBS7anegative_iamp21_path_pmc, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), "_.+")
  sample_name <- metadata_biobank$Anonymous_ID[grep(sample_name, metadata_biobank$Tumor_ID)]
  SBS7anegative_iamp21s_names_pmc <- c(SBS7anegative_iamp21s_names_pmc, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7anegative_iamp21s_pmc_gr <- sample_gr
  } else {
    SBS7anegative_iamp21s_pmc_gr <- c(SBS7anegative_iamp21s_pmc_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#Plot per sample ---------------------------------------------------------------
#Region 1
#Set parameters
for (sample in unique(SBS7apositive_iamp21s_stjude_gr$sample)){
  sample_gr <- SBS7apositive_iamp21s_stjude_gr[SBS7apositive_iamp21s_stjude_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[20e6:25e6])
  sample_df <- data.frame(sample = sample,
                          dataset = "StJude",
                          SBS7aStatus = "SBS7aPositive",
                          mean_coverage = mean_coverage,
                          region = "1")
  if (sample == unique(SBS7apositive_iamp21s_stjude_gr$sample)[1]){
    coverage_df <- sample_df
  } else {
    coverage_df <- rbind(coverage_df, sample_df)
  }
}
for (sample in unique(SBS7anegative_iamp21s_stjude_gr$sample)){
  sample_gr <- SBS7anegative_iamp21s_stjude_gr[SBS7anegative_iamp21s_stjude_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[20e6:25e6])
  sample_df <- data.frame(sample = sample,
                          dataset = "StJude",
                          SBS7aStatus = "SBS7aNegative",
                          mean_coverage = mean_coverage,
                          region = "1")
  coverage_df <- rbind(coverage_df, sample_df)
}
for (sample in unique(SBS7apositive_iamp21s_pmc_gr$sample)){
  sample_gr <- SBS7apositive_iamp21s_pmc_gr[SBS7apositive_iamp21s_pmc_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[20e6:25e6])
  sample_df <- data.frame(sample = sample,
                          dataset = "PMC",
                          SBS7aStatus = "SBS7aPositive",
                          mean_coverage = mean_coverage,
                          region = "1")
  coverage_df <- rbind(coverage_df, sample_df)
}
for (sample in unique(SBS7anegative_iamp21s_pmc_gr$sample)){
  sample_gr <- SBS7anegative_iamp21s_pmc_gr[SBS7anegative_iamp21s_pmc_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[20e6:25e6])
  sample_df <- data.frame(sample = sample,
                          dataset = "PMC",
                          SBS7aStatus = "SBS7aNegative",
                          mean_coverage = mean_coverage,
                          region = "1")
  coverage_df <- rbind(coverage_df, sample_df)
}

#Region 2
#Set parameters
for (sample in unique(SBS7apositive_iamp21s_stjude_gr$sample)){
  sample_gr <- SBS7apositive_iamp21s_stjude_gr[SBS7apositive_iamp21s_stjude_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[40e6:46709983])
  sample_df <- data.frame(sample = sample,
                          dataset = "StJude",
                          SBS7aStatus = "SBS7aPositive",
                          mean_coverage = mean_coverage,
                          region = "2")
  coverage_df <- rbind(coverage_df, sample_df)
}
for (sample in unique(SBS7anegative_iamp21s_stjude_gr$sample)){
  sample_gr <- SBS7anegative_iamp21s_stjude_gr[SBS7anegative_iamp21s_stjude_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[40e6:46709983])
  sample_df <- data.frame(sample = sample,
                          dataset = "StJude",
                          SBS7aStatus = "SBS7aNegative",
                          mean_coverage = mean_coverage,
                          region = "2")
  coverage_df <- rbind(coverage_df, sample_df)
}
for (sample in unique(SBS7apositive_iamp21s_pmc_gr$sample)){
  sample_gr <- SBS7apositive_iamp21s_pmc_gr[SBS7apositive_iamp21s_pmc_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[40e6:46709983])
  sample_df <- data.frame(sample = sample,
                          dataset = "PMC",
                          SBS7aStatus = "SBS7aPositive",
                          mean_coverage = mean_coverage,
                          region = "2")
  coverage_df <- rbind(coverage_df, sample_df)
}
for (sample in unique(SBS7anegative_iamp21s_pmc_gr$sample)){
  sample_gr <- SBS7anegative_iamp21s_pmc_gr[SBS7anegative_iamp21s_pmc_gr$sample == sample]
  sample_weighted_coverage <- coverage(sample_gr, weight = sample_gr$COPY_NUMBER)[["chr21"]]
  mean_coverage <- mean(as.vector(sample_weighted_coverage)[40e6:46709983])
  sample_df <- data.frame(sample = sample,
                          dataset = "PMC",
                          SBS7aStatus = "SBS7aNegative",
                          mean_coverage = mean_coverage,
                          region = "2")
  coverage_df <- rbind(coverage_df, sample_df)
}
#Check mean coverage
ggplot(coverage_df, aes(x = SBS7aStatus, y = mean_coverage, fill = SBS7aStatus)) +
  geom_boxplot() +
  facet_wrap(~region) +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_light() +
  geom_signif(comparisons = list(c("SBS7aNegative", "SBS7aPositive")),
              map_signif_level = T,
              test = "wilcox.test")
ggsave("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/iAMP21/Region_Coverage_Combined.pdf",
       width = 10, height = 5, units = "in")
