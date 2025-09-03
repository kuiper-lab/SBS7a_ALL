#Load packages
library(tidyverse)
library(cowplot)
library(ggsignif)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load data
bam_read_counts <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/deepsequencing/bam_read_counts_combined_28032023.tsv")
variants_list <- read_tsv("~/Projects/hypermutated_ALL/ANALYSES/deepsequencing/variants_list.tsv")

#Combine data
combined_data <- left_join(bam_read_counts, variants_list, by = c("#Chromosome" = "chr",
                                                                  "Position" = "start"))

#Get AF in R1 and Dx
combined_data$sample_R1_AF <- NA
combined_data$sample_Dx_AF <- NA
for (i in 1:nrow(combined_data)){
  combined_data$sample_Dx_AF[i] <- combined_data[i, paste0(combined_data$alt, "_AF_", paste0(str_remove(combined_data$sample, "P_"), "Dx"))[i]]
  combined_data$sample_R1_AF[i] <- combined_data[i, paste0(combined_data$alt, "_AF_", paste0(str_remove(combined_data$sample, "P_"), "R1"))[i]]
}
combined_data$sample_R1_reads <- NA
combined_data$sample_Dx_reads <- NA
for (i in 1:nrow(combined_data)){
  combined_data$sample_Dx_reads[i] <- combined_data[i, paste0(combined_data$alt, "_Number_of_reads_", paste0(str_remove(combined_data$sample, "P_"), "Dx"))[i]]
  combined_data$sample_R1_reads[i] <- combined_data[i, paste0(combined_data$alt, "_Number_of_reads_", paste0(str_remove(combined_data$sample, "P_"), "R1"))[i]]
}

#Look at depth
total_reads_cols <- grep("Total_reads", colnames(combined_data), value = T)
total_reads_df <- combined_data[, grep("Undetermined", total_reads_cols, invert = T, value = T)]
total_reads_df_long <- pivot_longer(total_reads_df, cols = 1:6, names_to = "sample", values_to = "total_reads")
total_reads_df_long$sample <- str_remove(total_reads_df_long$sample, "Total_reads_")
ggplot(total_reads_df_long) +
  geom_boxplot(aes(y = total_reads, x = sample, fill = sample)) +
  scale_fill_manual(values = cbp2) +
  theme_light()

#P0608 only
P0608_variants <- combined_data[combined_data$sample == "P0608",]

for (i in 1:nrow(P0608_variants)){
  alt_nuc <- P0608_variants$alt[i]
  alt_read_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_Number_of_reads"), colnames(P0608_variants), value = T), invert = T, value = T)
  new_cols <- paste0("alt_reads_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  alt_AF_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_AF"), colnames(P0608_variants), value = T), invert = T, value = T)
  new_AF_cols <- paste0("alt_AF_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  if (i == 1){
    P0608_variants[, new_cols] <- NA
  }
  P0608_variants[i, new_cols] <- P0608_variants[i, alt_read_cols]
  P0608_variants[i, new_AF_cols] <- P0608_variants[i, alt_AF_cols]
}

P0608_only_alt_AF <- P0608_variants[, c(grep("alt_AF", colnames(P0608_variants), value = T), "variant")]
P0608_only_alt_AF_long <- pivot_longer(P0608_only_alt_AF, cols = 1:6, names_to = "sample", values_to = "VAF")

#P0557 only
P0557_variants <- combined_data[combined_data$sample == "P0557",]

for (i in 1:nrow(P0557_variants)){
  alt_nuc <- P0557_variants$alt[i]
  alt_read_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_Number_of_reads"), colnames(P0557_variants), value = T), invert = T, value = T)
  new_cols <- paste0("alt_reads_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  alt_AF_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_AF"), colnames(P0557_variants), value = T), invert = T, value = T)
  new_AF_cols <- paste0("alt_AF_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  if (i == 1){
    P0557_variants[, new_cols] <- NA
  }
  P0557_variants[i, new_cols] <- P0557_variants[i, alt_read_cols]
  P0557_variants[i, new_AF_cols] <- P0557_variants[i, alt_AF_cols]
}

P0557_only_alt_AF <- P0557_variants[, c(grep("alt_AF", colnames(P0557_variants), value = T), "variant")]
P0557_only_alt_AF_long <- pivot_longer(P0557_only_alt_AF, cols = 1:6, names_to = "sample", values_to = "VAF")

#P0611 only
P0611_variants <- combined_data[combined_data$sample == "P0611",]

for (i in 1:nrow(P0611_variants)){
  alt_nuc <- P0611_variants$alt[i]
  alt_read_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_Number_of_reads"), colnames(P0611_variants), value = T), invert = T, value = T)
  new_cols <- paste0("alt_reads_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  alt_AF_cols <- grep("Undetermined",grep(paste0(alt_nuc, "_AF"), colnames(P0611_variants), value = T), invert = T, value = T)
  new_AF_cols <- paste0("alt_AF_", str_remove(alt_read_cols,paste0(alt_nuc, "_Number_of_reads_")))
  if (i == 1){
    P0611_variants[, new_cols] <- NA
  }
  P0611_variants[i, new_cols] <- P0611_variants[i, alt_read_cols]
  P0611_variants[i, new_AF_cols] <- P0611_variants[i, alt_AF_cols]
}

P0611_only_alt_AF <- P0611_variants[, c(grep("alt_AF", colnames(P0611_variants), value = T), "variant")]
P0611_only_alt_AF_long <- pivot_longer(P0611_only_alt_AF, cols = 1:6, names_to = "sample", values_to = "VAF")

#Make final plots
##P0608
P0608_only_alt_AF_long_simplified <- P0608_only_alt_AF_long
P0608_only_alt_AF_long_simplified$sample[P0608_only_alt_AF_long_simplified$sample == "alt_AF_P0608Dx"] <- "Dx"
P0608_only_alt_AF_long_simplified$sample[P0608_only_alt_AF_long_simplified$sample == "alt_AF_P0608R1"] <- "R"
P0608_only_alt_AF_long_simplified$sample[grep("alt", P0608_only_alt_AF_long_simplified$sample)] <- "C"
plot_P0608 <- ggplot(P0608_only_alt_AF_long_simplified, aes(x = sample, y = VAF)) +
  geom_jitter(height = 0, aes(col = sample)) +
  theme_bw() +
  ylim(0, 1) +
  scale_color_manual(values = cbp2) +
  ggtitle("Selection of 19 mutations") +
  geom_signif(comparisons = list(
    c("C", "Dx"),
    c("Dx", "R")), 
    map_signif_level=TRUE, test = "wilcox.test",
    y_position = c(0.8, 0.98),
    margin_top = 0)

#P0557
P0557_only_alt_AF_long_simplified <- P0557_only_alt_AF_long
P0557_only_alt_AF_long_simplified$sample[P0557_only_alt_AF_long_simplified$sample == "alt_AF_P0557Dx"] <- "Dx"
P0557_only_alt_AF_long_simplified$sample[P0557_only_alt_AF_long_simplified$sample == "alt_AF_P0557R1"] <- "R"
P0557_only_alt_AF_long_simplified$sample[grep("alt", P0557_only_alt_AF_long_simplified$sample)] <- "C"
plot_P0557 <- ggplot(P0557_only_alt_AF_long_simplified, aes(x = sample, y = VAF)) +
  geom_jitter(height = 0, aes(col = sample)) +
  theme_bw() +
  ylim(0, 1) +
  scale_color_manual(values = cbp2) +
  ggtitle("Selection of 15 mutations") +
  geom_signif(comparisons = list(
    c("C", "Dx"),
    c("Dx", "R")), 
    map_signif_level=TRUE, test = "wilcox.test",
    y_position = c(0.8, 0.98),
    margin_top = 0)

#P0611
P0611_only_alt_AF_long_simplified <- P0611_only_alt_AF_long
P0611_only_alt_AF_long_simplified$sample[P0611_only_alt_AF_long_simplified$sample == "alt_AF_P0611Dx"] <- "Dx"
P0611_only_alt_AF_long_simplified$sample[P0611_only_alt_AF_long_simplified$sample == "alt_AF_P0611R1"] <- "R"
P0611_only_alt_AF_long_simplified$sample[grep("alt", P0611_only_alt_AF_long_simplified$sample)] <- "C"
plot_P0611 <- ggplot(P0611_only_alt_AF_long_simplified, aes(x = sample, y = VAF)) +
  geom_jitter(height = 0, aes(col = sample)) +
  theme_bw() +
  ylim(0, 1) +
  scale_color_manual(values = cbp2) +
  ggtitle("Selection of 16 mutations") +
  geom_signif(comparisons = list(
    c("C", "Dx"),
    c("Dx", "R")), 
    map_signif_level=TRUE, test = "wilcox.test",
    y_position = c(0.8, 0.98),
    margin_top = 0)

#Save plot
pdf("~/Projects/hypermutated_ALL/RESULTS/plots/paper/DeepSequencing_vertical_significance.pdf", height = 10, width = 4)
plot_grid(plot_P0557, plot_P0608, plot_P0611, nrow = 3)
dev.off()
save(object = plot_P0557, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0557_deepseq_significance.rdata")
save(object = plot_P0608, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0608_deepseq_significance.rdata")
save(object = plot_P0611, file = "~/Projects/hypermutated_ALL/ANALYSES/MutationalTiming/P0611_deepseq_significance.rdata")
