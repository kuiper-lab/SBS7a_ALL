#Load libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(PCAtools)
library(apeglm)
library(vsn)
library(EnhancedVolcano)
library(hexbin)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)

#Load metadata and count matrix
load("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis/hyperdiploid_countmatrix.rdata") #Can be changed to different dataset
load("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis/hyperdiploid_metadata.rdata") #Can be changed to different dataset

#Make DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ sex + SBS7aStatus)

#Pre-filtering de dataset
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

#Transformation
#This transforms the data so analyses will not be influenced only by genes with a high count
vsd <- vst(dds, blind = FALSE)
meanSdPlot(assay(vsd))

#Run DE analysis
dds$SBS7aStatus <- relevel(dds$SBS7aStatus, ref = "SBS7aNegative")
dds <- DESeq(dds)
resultsNames(dds)
res05 <- results(dds, alpha=0.05, name="SBS7aStatus_SBS7aPositive_vs_SBS7aNegative")
summary(res05)
table(res05$padj < 0.05)
topGene <- rownames(res05)[which.min(res05$padj)]
plotCounts(dds, gene = topGene, intgroup="SBS7aStatus")

#Make volcano plot
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis/VolcanoPlot.pdf")
EnhancedVolcano(res05, lab = rownames(res05), x = 'log2FoldChange', y = "pvalue",
                pCutoffCol = "padj",
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0,10))
dev.off()

#Save list of differentially expressed genes
DE_genes <- res05[res05$padj < 0.05 & !is.na(res05$padj),]
DEgenes_df <- data.frame(gene = rownames(DE_genes),
                         direction = ifelse(DE_genes$log2FoldChange > 0, "up", "down"),
                         log2FC = DE_genes$log2FoldChange)
write.table(DEgenes_df, "~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis/DEgenes.tsv")


