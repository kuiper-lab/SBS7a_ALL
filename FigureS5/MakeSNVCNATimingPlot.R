#Load packages
library(tidyverse)
Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")

#Load data
germline_snvs <- read.table("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/GermlineGainedSNVs_UV_010.tsv")
somatic_snvs <- read.table("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/CNA_SNV_timing/SomaticGainedSNVs_UV_010.tsv")

#Combine data in 1 file
germline_snvs$`SNV Timing` <- "Germline SNVs"
somatic_snvs$`SNV Timing` <- "Somatic Mutations"
snvs_combined <- rbind(germline_snvs, somatic_snvs)

#Change names
snvs_combined$subtype[grep("hyperdiploid", snvs_combined$subtype, ignore.case = T)] <- "High Hyperdiploid"

#Remove samples with subtypes which are not informative
snvs_combined_selected <- snvs_combined[snvs_combined$subtype %in% c("iAMP21",
                                                                     "High Hyperdiploid"), ]

#Plot
CNA_SNV_timing_plot <- ggplot(snvs_combined_selected, aes(y = VAF, x = `SNV Timing`, fill = subtype)) +
  geom_violin(linewidth = 0.5) +
  scale_fill_manual(values = c(Tol_light[c(1,3)]), name = "Subtype") +
  facet_wrap(~subtype) +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "lightgray"),
        panel.border = element_rect(color = "lightgray", fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 15)) +
  ylab("VAF in tumor") +
  xlab("Chromosome 21 variants in gained regions")
ggsave(plot = CNA_SNV_timing_plot,
       filename = "~/surfdrive/Shared/Kuiper group/Relapsed_ALL/PROJECTS/UV SBS7a/Figures/Supplementary/CNA_SNV_timing.pdf",
       width = 10, height = 6, units = "in")
