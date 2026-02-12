#Load packages
library(tidyverse)
library(clusterProfiler)

#Load dataframes
iamp_df <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/GSEAiAMP21vsFusiondriven.tsv")
hyperdiploid_df <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/GSEAHyperdiploidvsFusiondriven.tsv")

#Look for overlap
intersect_df <- hyperdiploid_df[hyperdiploid_df$ID %in% iamp_df$ID,]
write.table(intersect_df, "~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/hyperdiploid_intersect.tsv")

#Compare genes
iamp_genes_df <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_06022026/DEgenesiAMP21vsFusiondriven.tsv")
hyperdiploid_genes_df <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_06022026/DEgenesHyperdiploidvsFusiondriven.tsv")
up_genes <- iamp_genes_df$gene[iamp_genes_df$log2FC >= 1 &
                                 iamp_genes_df$gene %in% hyperdiploid_genes_df$gene[hyperdiploid_genes_df$log2FC >= 1]]
up_genes <- str_remove(up_genes, ".+\\(") %>%
  str_remove("\\..+")
down_genes <- iamp_genes_df$gene[iamp_genes_df$log2FC <= -1 &
                                 iamp_genes_df$gene %in% hyperdiploid_genes_df$gene[hyperdiploid_genes_df$log2FC <= -1]]
down_genes <- str_remove(down_genes, ".+\\(") %>%
  str_remove("\\..+")

# Conduct over-representation analysis
ora_up <- enrichGO(gene = up_genes,
                   keyType = "ENSEMBL",
                   ont = "ALL",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   OrgDb = "org.Hs.eg.db", 
                   minGSSize = 3, 
                   maxGSSize = 1000)
ora_up_df <- as.data.frame(ora_up@result)
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/ora_upregulated_dotplot.pdf", width = 12)
dotplot(ora_up,
        x = "GeneRatio",
        showCategory = 50, # Choose how many categories should be shown
        font.size = 8, 
        title = "Top 50 enriched GO terms in upregulated genes") + 
  scale_y_discrete(labels = function(ora_up) str_wrap(ora_up, width = 100)) +
  theme(strip.text.y = element_text(size = 13),
        plot.title = element_text(face = "bold", size = 14))
dev.off()

ora_down <- enrichGO(gene = down_genes,
                     keyType = "ENSEMBL",
                     ont = "ALL",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "fdr",
                     OrgDb = "org.Hs.eg.db",
                     minGSSize = 3, 
                     maxGSSize = 1000)
ora_down_df <- as.data.frame(ora_down@result)
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/ora_downregulated_dotplot.pdf", width = 12)
dotplot(ora_down,
        x = "GeneRatio",
        showCategory = 50, # Choose how many categories should be shown
        font.size = 8, 
        title = "Top 50 enriched GO terms in downregulated genes") + 
  scale_y_discrete(labels = function(ora_up) str_wrap(ora_up, width = 100)) +
  theme(strip.text.y = element_text(size = 13),
        plot.title = element_text(face = "bold", size = 14))
dev.off()

#Save tables
ora_up_df$direction <- "up"
ora_down_df$direction <- "down"
ora_full_df <- rbind(ora_up_df,
                     ora_down_df)
write_tsv(ora_full_df, "~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/ORA_results.tsv")
