#Load packages
library(tidyverse)
library(ggsignif)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load mutational matrix
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_mut_mat.rdata")

#Get SBS7a positive samples
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/SBS7aPositiveSamples_SBS7aNoClusters.rdata")
SBS7_samples_ALL <- SBS7_samples_complete
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_SBS7_samples.rdata")

#Select only initial diagnoses
SBS7_samples_ALL_Dx <- grep("Dx",
                            SBS7_samples_ALL,
                            value = T)
SBS7_samples_ALCL_skincancer_Dx <- grep(pattern = "P0", SBS7_samples_complete, value = T, invert = T) %>%
  grep(pattern = "TI", value = T, invert = T)

#Get mutational loads
mut_loads <- as.data.frame(colSums(mut_mat))
colnames(mut_loads) <- "mutational_load"
mut_loads$dataset <- NA
mut_loads$dataset[rownames(mut_loads) %in% SBS7_samples_ALL_Dx] <- "SBS7a-positive BCP-ALL"
mut_loads$dataset[grep("P0\\d+$", rownames(mut_loads))] <- "SBS7a-negative BCP-ALL"
mut_loads$dataset[grepl("^PM", rownames(mut_loads)) &
                    rownames(mut_loads) %in% SBS7_samples_ALCL_skincancer_Dx] <- "SBS7a-positive ALCL"
mut_loads$dataset[rownames(mut_loads) %in%
                    SBS7_samples_ALCL_skincancer_Dx[grep("^PM",
                                                         SBS7_samples_ALCL_skincancer_Dx,
                                                         invert = T)]] <- "SBS7a-positive Skin Cancer"

#Select only mutational loads that will be used in this analysis
mut_loads_selected <- mut_loads[!is.na(mut_loads$dataset),]

#Set order
mut_loads_selected$dataset <- factor(mut_loads_selected$dataset, levels = c("SBS7a-negative BCP-ALL",
                                                                            "SBS7a-positive BCP-ALL",
                                                                            "SBS7a-positive ALCL",
                                                                            "SBS7a-positive Skin Cancer"))

#Plot
burden_plot <- ggplot(mut_loads_selected, aes(x = dataset, y = mutational_load, fill = dataset)) +
  geom_boxplot() +
  scale_y_log10(limits = c(200,4000000)) +
  geom_signif(comparisons = list(c("SBS7a-negative BCP-ALL", "SBS7a-positive BCP-ALL"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive Skin Cancer"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive ALCL")), 
              map_signif_level=TRUE, test = "wilcox.test",
              y_position = c(4.5, 6.3, 5.4), vjust = 0) +
  theme_light() +
  scale_fill_manual(values=c("lightgrey", cbp2[6], cbp2[1], cbp2[4]), name = "Tumor Type") +
  scale_x_discrete(
    element_blank(),
    labels = element_blank()
  ) +
  ylab(label = "Mutational Load") +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0))
save(object = burden_plot, file = "~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/mutationalloadplot.rdata")

#Calculate ratios
mut_loads_selected$mutational_load[mut_loads_selected$dataset == "SBS7a-negative BCP-ALL"] %>% mean() -> mean_ALL_no_UV
mut_loads_selected$mutational_load[mut_loads_selected$dataset == "SBS7a-positive BCP-ALL"] %>% mean() -> mean_ALL_UV
mut_loads_selected$mutational_load[mut_loads_selected$dataset == "SBS7a-positive Skin Cancer"] %>% mean() -> mean_Skin
mut_loads_selected$mutational_load[mut_loads_selected$dataset == "SBS7a-positive ALCL"] %>% mean() -> mean_ALCL
mean_ALL_UV / mean_ALL_no_UV
mean_Skin / mean_ALL_UV
