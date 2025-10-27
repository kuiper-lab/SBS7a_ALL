#Load packages
library(tidyverse)
library(ggsignif)
library(readxl)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load metadata
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/denovoExtraction/SBS7aPositiveSamples_SBS7aNoClusters.rdata")
SBS7_samples_ALL <- SBS7_samples_complete
load("~/Projects/hypermutated_ALL/ANALYSES/tumormutationburden/ALL_ALCL_skincancer_SBS7_samples.rdata")
metadata <- as.data.frame(read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx"))
metadata_subtype <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/metadata/B-ALL_subtypes_manual.xlsx")

#Get Transcriptional Strand Bias
load("~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/ALLstrandbiasNoClusters_multipletimepoints_nogermlinefilter_nosoftclippedbases_nosoftclippedbases_015_25042025.rdata")
SBS7aPositiveALL_TSB <- combined_strand_bias
recode_vector <- setNames(metadata$Coded_num, metadata$SKION)[unique(str_extract(SBS7aPositiveALL_TSB$group, "\\d\\d\\d+"))]
SBS7aPositiveALL_TSB$group <- str_replace_all(SBS7aPositiveALL_TSB$group, recode_vector)
SBS7a_samples_ALL_unique <- sort(SBS7_samples_ALL)[!duplicated(unlist(str_extract_all(sort(SBS7_samples_ALL), "P\\d\\d\\d\\d+")))]
SBS7aPositiveALL_TSB <- SBS7aPositiveALL_TSB[SBS7aPositiveALL_TSB$group %in% SBS7a_samples_ALL_unique,]

load("~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/SBS7aNegativeBCPALLstrandbias_18042025.rdata")
SBS7aNegativeALL_TSB <- SBS7aNegativeALL_TSB[SBS7aNegativeALL_TSB$group %in% metadata_subtype$Anonymous,]
SBS7aNegativeALL_TSB <- SBS7aNegativeALL_TSB[!SBS7aNegativeALL_TSB$group %in% str_extract(SBS7_samples_ALL, "P\\d+"),]

load("~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/ALCLstrandbias_14082024.rdata")
SBS7aPositiveALCL_TSB <- ALCL_strand_bias[ALCL_strand_bias$group %in% SBS7_samples_complete,]

SBS7aPositiveSkinCancer_TSB <- read.csv("~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/SkinCancerTSB_SBS7aPositive.csv",
                                        row.names = 1)

#Combine dataframes
SBS7aPositiveALL_TSB$TumorType <- "SBS7a-positive BCP-ALL"
SBS7aNegativeALL_TSB$TumorType <- "SBS7a-negative BCP-ALL"
SBS7aPositiveALCL_TSB$TumorType <- "SBS7a-positive ALCL"
SBS7aPositiveSkinCancer_TSB$TumorType <- "SBS7a-positive Skin Cancer"
combined_TSB <- rbind(SBS7aPositiveALL_TSB,
                      SBS7aNegativeALL_TSB,
                      SBS7aPositiveALCL_TSB,
                      SBS7aPositiveSkinCancer_TSB)
combined_TSB$TumorType <- factor(combined_TSB$TumorType, levels = c("SBS7a-negative BCP-ALL",
                                                                    "SBS7a-positive BCP-ALL",
                                                                    "SBS7a-positive ALCL",
                                                                    "SBS7a-positive Skin Cancer"))

#Get C>T mutations
combined_CT <- combined_TSB[combined_TSB$type == "C>T",]
TSB_plot <- ggplot(combined_CT, aes(x = TumorType, y = ratio, fill = TumorType)) +
  geom_boxplot() +
  labs(y = "Transcribed / Untranscribed") +
  theme_light() +
  geom_signif(comparisons = list(c("SBS7a-positive BCP-ALL", "SBS7a-positive Skin Cancer"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive ALCL"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-negative BCP-ALL")
  ), 
  map_signif_level=TRUE, test = "wilcox.test",
  y_position = c(1.3, 1.6, 1.9)) +
  scale_fill_manual(values=c("lightgrey", cbp2[6], cbp2[1], cbp2[4]), name = "Tumor Type") +
  scale_x_discrete(
    element_blank(),
    labels = element_blank()
  ) +
  ylim(c(0.4, 2)) +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0))
save(object = TSB_plot, file = "~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/TSB_plot.rdata")
TSB_plot_logtransformed <- ggplot(combined_CT, aes(x = TumorType, y = log2(ratio), fill = TumorType)) +
  geom_boxplot() +
  labs(y = expression(log[2](Transcribed / Untranscribed))) +
  theme_light() +
  geom_signif(comparisons = list(c("SBS7a-positive BCP-ALL", "SBS7a-positive Skin Cancer"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-positive ALCL"),
                                 c("SBS7a-positive BCP-ALL", "SBS7a-negative BCP-ALL")
  ), 
  map_signif_level=TRUE, test = "wilcox.test",
  y_position = c(0.3, 0.6, 0.9)) +
  scale_fill_manual(values=c("lightgrey", cbp2[6], cbp2[1], cbp2[4]), name = "Tumor Type") +
  scale_x_discrete(
    element_blank(),
    labels = element_blank()
  ) +
  theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0))
save(object = TSB_plot_logtransformed, file = "~/Projects/hypermutated_ALL/ANALYSES/melanoma_comparison/TranscriptionalStrandBias/TSB_logtransformed_plot.rdata")

#Perform statistical analysis
wilcox.test(log2(combined_CT$ratio[combined_CT$TumorType == "SBS7a-positive BCP-ALL"]),
            log2(combined_CT$ratio[combined_CT$TumorType == "SBS7a-negative BCP-ALL"]))
wilcox.test(log2(combined_CT$ratio[combined_CT$TumorType == "SBS7a-positive BCP-ALL"]),
            log2(combined_CT$ratio[combined_CT$TumorType == "SBS7a-positive Skin Cancer"]))
