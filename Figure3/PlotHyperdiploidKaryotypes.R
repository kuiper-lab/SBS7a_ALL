#Load packages -----------------------------------------------------------------
library(tidyverse)
library(pheatmap)
library(readxl)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load data----------------------------------------------------------------------
#Load CNV data
pmc_data <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/Hyperdiploid/Karyotype/hyperdiploid_overview_PMC_all_chromosomes.tsv")
stjude_data <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/Hyperdiploid/Karyotype/hyperdiploid_overview_StJude_all_chromosomes.tsv")

#Get SBS7a annotations
SBS7aPositive_pmc <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aPositiveHyperdiploids_karyotypes.xlsx")
SBS7aNegative_pmc <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aNegativeHyperdiploids_karyotypes.xlsx")
SBS7aNegative_StJude <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/StJude/Hyperdiploid_SBS7aNegativeComplete.txt",
                                 col_names = F)
SBS7aPositive_StJude <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/StJude/Hyperdiploid_SBS7aPositiveComplete.txt",
                                 col_names = F)

#Remove some samples
pmc_data_clean <- pmc_data[str_remove(pmc_data$Sample, "_.+") %in%
                             SBS7aNegative_pmc$Tumor_ID |
                             pmc_data$Sample %in% SBS7aPositive_pmc$Tumor_ID,]
colnames(pmc_data_clean) <- str_remove(colnames(pmc_data_clean), "_copynumber")
write_tsv(pmc_data_clean, "~/surfdrive/Shared/Kuiper group/Relapsed_ALL/ANALYSES/SBS7a/DATA/CNV/PMC/hyperdiploid_overview_PMC_all_chromosomes_filtered.tsv")
stjude_data_clean <- stjude_data[stjude_data$Sample %in%
                                   c(SBS7aNegative_StJude$X1,
                                     SBS7aPositive_StJude$X1), ]
colnames(stjude_data_clean) <- str_remove(colnames(stjude_data_clean), "_copynumber")
write_tsv(stjude_data_clean, "~/surfdrive/Shared/Kuiper group/Relapsed_ALL/ANALYSES/SBS7a/DATA/CNV/StJude/hyperdiploid_overview_StJude_all_chromosomes_filtered.tsv")

#Get full SBS7a annotations
full_positive_list <- c(SBS7aPositive_pmc$Tumor_ID, SBS7aPositive_StJude$X1)
full_negative_list <- c(SBS7aNegative_pmc$Tumor_ID, SBS7aNegative_StJude$X1)

#Make binary karyotypes --------------------------------------------------------
chromosomes <- paste0("chr", 1:22)
coded_karyotypes_pmc <- pmc_data_clean[, chromosomes]
coded_karyotypes_pmc <- apply(coded_karyotypes_pmc, 2, as.numeric)
rownames(coded_karyotypes_pmc) <- pmc_data_clean$Sample
coded_karyotypes_pmc[!is.na(coded_karyotypes_pmc)] <- 1
coded_karyotypes_pmc[is.na(coded_karyotypes_pmc)] <- 0

coded_karyotypes_stjude <- stjude_data_clean[, chromosomes]
coded_karyotypes_stjude <- apply(coded_karyotypes_stjude, 2, as.numeric)
rownames(coded_karyotypes_stjude) <- stjude_data_clean$Sample
coded_karyotypes_stjude[!is.na(coded_karyotypes_stjude)] <- 1
coded_karyotypes_stjude[is.na(coded_karyotypes_stjude)] <- 0

#Combine datasets
coded_karyotypes <- rbind(coded_karyotypes_pmc, coded_karyotypes_stjude)

#Make heatmap of chromosomes ---------------------------------------------------
mat <- coded_karyotypes[, paste0("chr", 1:22)]
anno <- data.frame("SBS7aStatus" = rownames(coded_karyotypes) %in% full_positive_list,
                   "Origin" = grepl("SJ", rownames(coded_karyotypes)))
anno$SBS7aStatus[anno$SBS7aStatus] <- "SBS7aPositive"
anno$SBS7aStatus[anno$SBS7aStatus == "FALSE"] <- "SBS7aNegative"
anno$Origin[anno$Origin] <- "StJude"
anno$Origin[anno$Origin == "FALSE"] <- "Biobank"
rownames(anno) <- rownames(mat)
dev.off()
pdf("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/ANALYSES/SBS7a/RESULTS/Figures/CNAcomparison/KaryotypeHyperdiploidHeatmap_manual.pdf",
    width = 10)
pheatmap(mat,
         annotation_row = anno,
         cluster_cols = F,
         show_rownames = F,
         color = c("lightgrey", "firebrick"),
         border_color = NA,
         annotation_colors = list("Origin" = c("Biobank" = "#F5793A",
                                               "StJude" = "#85C0F9"),
                                  "SBS7aStatus" = c("SBS7aPositive" = cbp2[4],
                                                    "SBS7aNegative" = "lightgrey")))
dev.off()

#Get Moorman prognoses ---------------------------------------------------------
moorman_prognosis <- c(rep(NA, nrow(mat)))
moorman_prognosis[rowSums(mat[, c("chr17", "chr18")]) == 2] <- "Good Risk"
moorman_prognosis[rowSums(mat[, c("chr17", "chr18")]) == 1 &
                    rowSums(mat[, c("chr5", "chr20")]) >= 1] <- "Poor Risk"
moorman_prognosis[rowSums(mat[, c("chr17", "chr18")]) == 1 &
                    rowSums(mat[, c("chr5", "chr20")]) == 0] <- "Good Risk"
moorman_prognosis[rowSums(mat[, c("chr17", "chr18")]) == 0] <- "Poor Risk"
anno$Moorman <- moorman_prognosis

sum(anno$Moorman == "Good Risk" & anno$SBS7aStatus == "SBS7aPositive")/
  sum(anno$SBS7aStatus == "SBS7aPositive")
sum(anno$Moorman == "Good Risk" & anno$SBS7aStatus == "SBS7aNegative")/
  sum(anno$SBS7aStatus == "SBS7aNegative")
ggplot(anno) +
  geom_bar(aes(x = SBS7aStatus, fill = Moorman), position = "fill") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)) +
  scale_x_discrete(labels=c("SBS7aNegative" = "SBS7a-negative",
                            "SBS7aPositive" = "SBS7a-positive"))
ggsave("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/ANALYSES/SBS7a/RESULTS/Figures/CNAcomparison/MoormanPrognosis.pdf",
       width = 3, height = 4, units = "in")
fisher.test(anno$SBS7aStatus, anno$Moorman)

#Plot percentage gained per chromosome -----------------------------------------
#Get percentage gained per SBS7a status
mat_SBS7aPositive <- mat[anno$SBS7aStatus == "SBS7aPositive",]
mat_SBS7aNegative <- mat[anno$SBS7aStatus == "SBS7aNegative",]
gain_percentage_df <- data.frame(SBS7aPositive = 100 * colSums(mat_SBS7aPositive) / nrow(mat_SBS7aPositive),
                                 SBS7aNegative = 100 * colSums(mat_SBS7aNegative) / nrow(mat_SBS7aNegative),
                                 Chromosome = colnames(mat)
                                 )
gain_percentage_df_long <- pivot_longer(gain_percentage_df,
                                        c(SBS7aPositive,SBS7aNegative),
                                        names_to = "SBS7aStatus",
                                        values_to = "Percentage Gained")
gain_percentage_df_long$Chromosome <- factor(gain_percentage_df_long$Chromosome,
                                             levels = colnames(mat))

#Plot
ggplot(gain_percentage_df_long, aes(x = Chromosome,
                                    y = `Percentage Gained`,
                                    fill = SBS7aStatus)) +
  geom_col(position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = c("lightgrey", cbp2[4]))
ggsave("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/ANALYSES/SBS7a/RESULTS/Figures/CNAcomparison/KaryotypeHyperdiploidGainPercentage.pdf",
    width = 10, height = 2, units = "in")


