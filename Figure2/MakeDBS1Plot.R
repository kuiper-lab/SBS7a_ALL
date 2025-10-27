#Load packages
library(tidyverse)
library(ggsignif)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load snp mutational matrix
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_mut_mat.rdata")

#Load mnp mutatational matrix
load("~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/DBS1/ALL_ALCL_skincancer_mnp_mut_mat.rdata")

#Make CT/CC>TT dataframe
colnames(mut_mat)[!colnames(mut_mat) %in% colnames(mut_mat_mnp)]
colnames(mut_mat_mnp)[!colnames(mut_mat_mnp) %in% colnames(mut_mat)]
samples_intersect <- colnames(mut_mat)
mut_mat_CT <- mut_mat[grep("\\[C>T\\]", rownames(mut_mat)),]
CT_CCTT_df <- data.frame(sample = colnames(mut_mat),
                         CT = colSums(mut_mat_CT)[colnames(mut_mat)],
                         CCTT = colSums(mut_mat_mnp["CC_TT", colnames(mut_mat)]))

#Get SBS7a positive samples
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/SBS7aPositiveSamples_SBS7aNoClusters.rdata")
SBS7_samples_ALL <- SBS7_samples_complete
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_SBS7_samples.rdata")

#Select only initial diagnoses
SBS7a_samples_ALL_unique <- sort(SBS7_samples_ALL)[!duplicated(unlist(str_extract_all(sort(SBS7_samples_ALL), "P\\d\\d\\d\\d+")))]
SBS7_samples_ALCL_skincancer_Dx <- grep(pattern = "P0", SBS7_samples_complete, value = T, invert = T) %>%
  grep(pattern = "TI", value = T, invert = T)

#Get mutational loads
CT_CCTT_df$dataset <- NA
CT_CCTT_df$dataset[rownames(CT_CCTT_df) %in% SBS7a_samples_ALL_unique] <- "SBS7a-positive BCP-ALL"
CT_CCTT_df$dataset[grep("P0\\d+$", rownames(CT_CCTT_df))] <- "SBS7a-negative BCP-ALL"
CT_CCTT_df$dataset[grepl("^PM", rownames(CT_CCTT_df)) &
                    rownames(CT_CCTT_df) %in% SBS7_samples_ALCL_skincancer_Dx] <- "SBS7a-positive ALCL"
CT_CCTT_df$dataset[rownames(CT_CCTT_df) %in%
                    SBS7_samples_ALCL_skincancer_Dx[grep("^PM",
                                                         SBS7_samples_ALCL_skincancer_Dx,
                                                         invert = T)]] <- "SBS7a-positive Skin Cancer"

#Select only samples that will be used in this analysis
CT_CCTT_df_selected <- CT_CCTT_df[!is.na(CT_CCTT_df$dataset),]

#Set order
CT_CCTT_df_selected$dataset <- factor(CT_CCTT_df_selected$dataset, levels = c("SBS7a-negative BCP-ALL",
                                                                              "SBS7a-positive BCP-ALL",
                                                                              "SBS7a-positive ALCL",
                                                                              "SBS7a-positive Skin Cancer"))

#Plot
DBS1_plot <- ggplot(CT_CCTT_df_selected,
                    aes(x = CT,
                        y = CCTT,
                        color = dataset)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_light() +
  labs(x="# C>T mutations",
       y="# CC>TT mutations") +
  scale_color_manual(values=c("lightgrey", cbp2[6], cbp2[1], cbp2[4]), name = "Tumor Type")
save(object = DBS1_plot, file = "~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/DBS1/CT_CCTT_plot.rdata")

#Calculate ratio
CT_CCTT_df_selected$ratio <- CT_CCTT_df_selected$CCTT/CT_CCTT_df_selected$CT
ggplot(CT_CCTT_df_selected, aes(x = dataset, y = ratio, fill = dataset)) +
  geom_boxplot() +
  labs(y = "CC>TT / C>T") +
  theme_light() +
  geom_signif(comparisons = list(c("SBS7a-negative BCP-ALL", "SBS7a-positive BCP-ALL"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive Skin Cancer"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive ALCL")
                                 ), 
              map_signif_level=TRUE, test = "wilcox.test") +
  scale_fill_manual(values=c("lightgrey", cbp2[6], cbp2[1], cbp2[4]), name = "Tumor Type") +
  ylim(c(0, 0.13))
wilcox.test(CT_CCTT_df_selected$ratio[CT_CCTT_df_selected$dataset == "SBS7a-positive BCP-ALL"],
            CT_CCTT_df_selected$ratio[CT_CCTT_df_selected$dataset == "SBS7a-positive Skin Cancer"])
wilcox.test(CT_CCTT_df_selected$ratio[CT_CCTT_df_selected$dataset == "SBS7a-positive BCP-ALL"],
            CT_CCTT_df_selected$ratio[CT_CCTT_df_selected$dataset == "SBS7a-negative BCP-ALL"])
CT_CCTT_df_selected %>%
  group_by(dataset) %>%
  summarise(avg=mean(ratio))

#Perform statistical test on number of CC>TT mutations
wilcox.test(CT_CCTT_df_selected$CCTT[CT_CCTT_df_selected$dataset == "SBS7a-positive BCP-ALL"],
            CT_CCTT_df_selected$CCTT[CT_CCTT_df_selected$dataset == "SBS7a-positive Skin Cancer"])
