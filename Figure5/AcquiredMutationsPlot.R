#Load packages
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(tidyverse)
library(readxl)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
library(Polychrome)

#Load absolute contributions
absolute_contributions <- read.csv("~/Projects/hypermutated_ALL/RESULTS/tables/AbsoluteContributions_SBS7aALLClusters_ReuseSignatures.csv", row.names = 1)
colnames(absolute_contributions) <- str_replace(colnames(absolute_contributions), "X", "P_")

#Get clusters
absolute_contributions_clusters <- absolute_contributions[, grep("cluster", colnames(absolute_contributions))]

ClusterCountFirstAppearance <- data.frame(patient = c("P0608", "P0608", "P0608", "P0608",
                                                      "P0611", "P0611",
                                                      "P0557", "P0557"),
                                          Tumor = c("Dx", "R1", "R2", "R3",
                                                    "Dx", "R1",
                                                    "Dx", "R1"),
                                          Anonym = c("P0608", "P0608", "P0608", "P0608",
                                                     "P0611", "P0611",
                                                     "P0557", "P0557"))

absolute_contributions_clusters[, grep("P0608", colnames(absolute_contributions_clusters))]
ClusterCountFirstAppearance$SBS7a  <- 0
ClusterCountFirstAppearance$burden  <- 0

#Fill P0608
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0608" &
                              ClusterCountFirstAppearance$Tumor == "Dx"] <- absolute_contributions_clusters["SBS7a", "P0608_cluster12b"] +
  absolute_contributions_clusters["SBS7a", "P0608_cluster5"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0608" &
                                    ClusterCountFirstAppearance$Tumor == "Dx"] <- colSums(absolute_contributions_clusters)["P0608_cluster12b"] +
  colSums(absolute_contributions_clusters)["P0608_cluster5"]
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0608" &
                                    ClusterCountFirstAppearance$Tumor == "R1"] <- absolute_contributions_clusters["SBS7a", "P0608_cluster4"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0608" &
                                     ClusterCountFirstAppearance$Tumor == "R1"] <- colSums(absolute_contributions_clusters)["P0608_cluster4"]
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0608" &
                                    ClusterCountFirstAppearance$Tumor == "R2"] <- absolute_contributions_clusters["SBS7a", "P0608_cluster2"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0608" &
                                     ClusterCountFirstAppearance$Tumor == "R2"] <- apply(absolute_contributions_clusters, 2, sum)["P0608_cluster2"]
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0608" &
                                    ClusterCountFirstAppearance$Tumor == "R3"] <- 0
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0608" &
                                     ClusterCountFirstAppearance$Tumor == "R3"] <- 0

#Fill P0611
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0611" &
                                    ClusterCountFirstAppearance$Tumor == "Dx"] <- absolute_contributions_clusters["SBS7a", "P0611_cluster2"] +
  absolute_contributions_clusters["SBS7a", "P0611_cluster3"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0611" &
                                     ClusterCountFirstAppearance$Tumor == "Dx"] <- apply(absolute_contributions_clusters, 2, sum)["P0611_cluster2"] +
  apply(absolute_contributions_clusters, 2, sum)["P0611_cluster3"]
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0611" &
                                    ClusterCountFirstAppearance$Tumor == "R1"] <- absolute_contributions_clusters["SBS7a", "P0611_cluster1"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0611" &
                                     ClusterCountFirstAppearance$Tumor == "R1"] <- apply(absolute_contributions_clusters, 2, sum)["P0611_cluster1"]

#Fill P0557
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0557" &
                                    ClusterCountFirstAppearance$Tumor == "Dx"] <- absolute_contributions_clusters["SBS7a", "P0557_cluster3"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0557" &
                                     ClusterCountFirstAppearance$Tumor == "Dx"] <- apply(absolute_contributions_clusters, 2, sum)["P0557_cluster3"]
ClusterCountFirstAppearance$SBS7a[ClusterCountFirstAppearance$patient == "P0557" &
                                    ClusterCountFirstAppearance$Tumor == "R1"] <- absolute_contributions_clusters["SBS7a", "P0557_cluster2"]
ClusterCountFirstAppearance$burden[ClusterCountFirstAppearance$patient == "P0557" &
                                     ClusterCountFirstAppearance$Tumor == "R1"] <- apply(absolute_contributions_clusters, 2, sum)["P0557_cluster2"]

#Make barplots per patient
ClusterCountFirstAppearance$Other <- ClusterCountFirstAppearance$burden - ClusterCountFirstAppearance$SBS7a
ClusterCountFirstAppearance_long <- pivot_longer(ClusterCountFirstAppearance,
                                                 c(SBS7a, Other),
                                                 names_to = "Signature",
                                                 values_to = "Relative Contribution of Acquired Mutations")
P0608_relative <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0608",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature), position = "fill") +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ggtitle("P0608")
P0557_relative <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0557",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature), position = "fill") +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ggtitle("P0557")
P0611_relative <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0611",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature), position = "fill") +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ggtitle("P0611")

#Absolute contributions
P0608_absolute <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0608",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature)) +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ylab("Absolute Contribution of Acquired Mutations") +
  ggtitle("P0608")
P0557_absolute <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0557",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature)) +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ylab("Absolute Contribution of Acquired Mutations") +
  ggtitle("P0557")
P0611_absolute <- ggplot(ClusterCountFirstAppearance_long[ClusterCountFirstAppearance_long$patient == "P0611",]) +
  geom_col(aes(x = Tumor, y = `Relative Contribution of Acquired Mutations`, fill = Signature)) +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme_classic() +
  ylab("Absolute Contribution of Acquired Mutations") +
  ggtitle("P0611")

#Save plots
save(object = P0608_absolute, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0608_absolute.rdata")
save(object = P0608_relative, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0608_relative.rdata")
save(object = P0557_absolute, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0557_absolute.rdata")
save(object = P0557_relative, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0557_relative.rdata")
save(object = P0611_absolute, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0611_absolute.rdata")
save(object = P0611_relative, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0611_relative.rdata")

