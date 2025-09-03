#Load packages
library(tidyverse)
library(ggvenn)
library(venneuler)
library(eulerr)

#Load differentially expressed genes
iAMP21_unstranded <- read.table("~/Projects/hypermutated_ALL/ANALYSES/StJude/RNAseq/iAMP21/DEgenesiAMP21Unstranded.tsv")
iAMP21_stranded <- read.table("~/Projects/hypermutated_ALL/ANALYSES/StJude/RNAseq/iAMP21/DEgenesiAMP21Stranded.tsv")
hyperdiploid_stjude <- read.table("~/Projects/hypermutated_ALL/ANALYSES/StJude/RNAseq/hyperdiploid/DEgenesHyperdiploid.tsv")
hyperdiploid_pmc <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/DEgenesHyperdiploid.tsv")
iAMP21_pmc <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/DEgenesiAMP21.tsv")

#Make venn diagram
venn_list <- list("iAMP21_unstranded" = as.character(iAMP21_unstranded$gene),
                  "iAMP21_stranded" = as.character(iAMP21_stranded$gene),
                  "hyperdiploid_StJude" = as.character(hyperdiploid_stjude$gene),
                  "hyperdiploid_PMC" = as.character(hyperdiploid_pmc$gene),
                  "iAMP21_PMC" = as.character(iAMP21_pmc$gene))
ggvenn(venn_list)

#Plot complete venn diagram
plot(venn(venn_list))

#Take direction into account and compare within subtypes
##Hyperdiploid
plot(venn(list("hyperdiploid_StJude" = as.character(hyperdiploid_stjude$gene),
               "hyperdiploid_PMC" = as.character(hyperdiploid_pmc$gene))))
hyperdiploid_pmc$gene[hyperdiploid_pmc$gene %in% hyperdiploid_stjude$gene]

##iAMP21
plot(venn(list("iAMP21_unstranded" = as.character(iAMP21_unstranded$gene),
               "iAMP21_stranded" = as.character(iAMP21_stranded$gene),
               "iAMP21_PMC" = as.character(iAMP21_pmc$gene))))
iamp_stranded_genes <- iAMP21_pmc$gene[iAMP21_pmc$gene %in% iAMP21_stranded$gene]
iamp_unstranded_genes <- iAMP21_pmc$gene[iAMP21_pmc$gene %in% iAMP21_unstranded$gene]
iAMP21_pmc[iAMP21_pmc$gene %in% iamp_stranded_genes,]
iAMP21_stranded[iAMP21_stranded$gene %in% iamp_stranded_genes,]
iAMP21_pmc[iAMP21_pmc$gene %in% iamp_unstranded_genes,]
iAMP21_unstranded[iAMP21_unstranded$gene %in% iamp_unstranded_genes,]

#Make final Venn diagrams
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/Hyperdiploid_Venn.pdf")
plot(venn(list("Hyperdiploid\nSt. Jude Children's\nResearch Hospital" = as.character(hyperdiploid_stjude$gene),
"High Hyperdiploid\nPrincess Màxima Center\nfor Pediatric Oncology" = as.character(hyperdiploid_pmc$gene))))
dev.off()
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/iAMP21_Venn.pdf")
plot(venn(list("iAMP21 Stranded\nSt. Jude Children's\nResearch Hospital" = paste0(as.character(iAMP21_stranded$gene), iAMP21_stranded$direction),
               "iAMP21 Unstranded\nSt. Jude Children's\nResearch Hospital" = paste0(as.character(iAMP21_unstranded$gene), iAMP21_unstranded$direction),
               "iAMP21\nPrincess Màxima Center\nfor Pediatric Oncology" = paste0(as.character(iAMP21_pmc$gene), iAMP21_pmc$direction))))
dev.off()

#Make table
iAMPdf_1 <- full_join(iAMP21_pmc, iAMP21_stranded,
                      by = "gene", suffix = c(" iAMP21 Princess Màxima Center for Pediatric Oncology",
                                              " iAMP21 Stranded St. Jude Children's Research Hospital"))
iAMP_df <- full_join(iAMPdf_1, iAMP21_unstranded, by = "gene")
colnames(iAMP_df)[c(6,7)] <- paste0(colnames(iAMP_df)[c(6,7)],
                                    " iAMP21 Unstranded St. Jude Children's Research Hospital")
