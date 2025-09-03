#Load packages
library(tidyverse)
library(cowplot)

#Load rdata
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0429.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0557_strict.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0621.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0628_strict.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0629.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/plots/P0633.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/DiagnosisMutationMatrices/SBS7aContribution_selected.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/DiagnosisMutationMatrices_strict/SBS7aContribution.rdata")

#Make absolute contribution plots
relativecontribution_plot_a <- ggplot(SBS7a_contribution_selected[SBS7a_contribution_selected$sample %in% c("P0429Dx",
                                                                                                            "P0629Dx",
                                                                                                            "P0633Dx"),],
                                         aes(x = vaf,
                                             y = contribution,
                                             group = sample,
                                             col = sample)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("") +
  ylab("Relative Contribution SBS7a") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 8))
relativecontribution_plot_b <- ggplot(rbind(SBS7a_contribution_selected[SBS7a_contribution_selected$sample %in% c("P0621Dx"),],
                                            SBS7a_contribution[SBS7a_contribution$sample %in% c("P0628Dx",
                                                                                                "P0557Dx"),]),
                                      aes(x = vaf,
                                          y = contribution,
                                          group = sample,
                                          col = sample)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  xlab("") +
  ylab("Relative Contribution SBS7a") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 8))

#Make complete plot
pdf("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/SubclonalMutationsMainFigureCategories_strict.pdf",
    width = 12)
plot_grid(fullplot_P0629, relativecontribution_plot_a,
          fullplot_P0621, relativecontribution_plot_b,
          ncol = 2, rel_widths = c(3,1), labels = c("A", "",
                                                    "B", ""))
dev.off()

#Make supplementary plot
sample_plots_sup <- plot_grid(fullplot_P0429,
                              fullplot_P0633,
                              fullplot_P0557,
                              fullplot_P0628,
                              ncol = 1)
pdf("~/Projects/hypermutated_ALL/ANALYSES/subclonal_mutations/SubclonalMutationsSupplementaryFigure_strict.pdf",
    width = 9)
sample_plots_sup
dev.off()
