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
library(readxl)

#Load metadata and count matrix
load("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/hyperdiploid_countmatrix.rdata")
load("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/hyperdiploid_metadata.rdata")

#Load NER genes
NER_genes <- scan("~/Projects/hypermutated_ALL/ANALYSES/NER/KEGG_NER_symbols.txt",
                  character(), quote = "", sep = "\n")

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

#Plot heatmap of NER genes
mat  <- assay(vsd)[ NER_genes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("SBS7aStatus","sex")])
dev.off()
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_05112024/HyperdiploidHeatmapNER.pdf")
pheatmap(mat, annotation_col = anno)
dev.off()
