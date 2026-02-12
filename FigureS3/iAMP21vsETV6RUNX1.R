#General---------------------------------------------------------------------
#Load libraries
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(readxl)

#Data preprocessing------------------------------------------------------------
#Load data
load("~/Projects/hypermutated_ALL/DATA/RNAseq/count_data_ALLsamples_PMC.Rda")

#Combine data
count_data_df_selection <- count_data_df[,grep("rawCounts", colnames(count_data_df))]
count_data_df_selection$Ensembl <- str_remove(unlist(lapply(str_split(rownames(count_data_df_selection), "\\("), '[', 2)), "\\)")
combined_data <- count_data_df_selection
colnames(combined_data)
combined_data_selection <- combined_data
all_patients <- str_remove(colnames(combined_data_selection)[grep("rawCounts", colnames(combined_data_selection))],
                           "_rawCounts")
colnames(combined_data_selection) <- str_remove(colnames(combined_data_selection), "_rawCounts")

#Make matrix
count_matrix <- as.matrix(combined_data_selection)
count_matrix_new <- count_matrix
count_matrix_numeric <- apply(count_matrix_new, 2, as.numeric)
rownames(count_matrix_numeric) <- rownames(count_matrix)
count_matrix <- count_matrix_numeric
count_matrix <- count_matrix[,!grepl("R", colnames(count_matrix))]

#Read metadata
load("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_07012026/Metadata.rdata")
coldata <- anno
count_matrix <- count_matrix[, rownames(coldata)]

#Remove high hyperdiploid B-ALL
metadata_subtypes <- read.delim("~/Projects/hypermutated_ALL/DATA/RNAseq/RNAseqSampleOverview.csv", sep = ",")
coldata <- coldata[!rownames(coldata) %in% metadata_subtypes$SKION[grep("hyperdiploid", metadata_subtypes$Sub.type, ignore.case = T)],]
count_matrix <- count_matrix[, rownames(coldata)]

#Make DEseq dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ Sex + Group)

#Pre-filtering de dataset
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

#Run DE analysis----------------------------------------------------------------
dds$Group <- relevel(dds$Group, ref = "Fusion-driven")
dds <- DESeq(dds)
resultsNames(dds)
res05 <- results(dds, alpha=0.05, name="Group_Aneuploid_vs_Fusion.driven")
res05$symbol <- str_remove(rownames(res05), " \\(.+")
res05$ENSEMBL <- str_remove(rownames(res05), ".+\\(") %>%
  str_remove("\\..+")

#Save list of differentially expressed genes
DE_genes <- res05[res05$padj < 0.05 & !is.na(res05$padj),]
DEgenes_df <- data.frame(gene = rownames(DE_genes),
                         direction = ifelse(DE_genes$log2FoldChange > 0, "up", "down"),
                         log2FC = DE_genes$log2FoldChange)
write.table(DEgenes_df, "~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_06022026/DEgenesiAMP21vsFusiondriven.tsv")

#Perform gene set enrichment analysis but get all GO terms ---------------------
# we want the log2 fold change 
original_gene_list <- res05$log2FoldChange
# name the vector
names(original_gene_list) <- res05$ENSEMBL
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# perform gse
gse <- gseGO(geneList=gene_list, 
                         ont ="ALL", 
                         keyType = "ENSEMBL", 
                         by = "fgsea",
                         #nPerm = 1000, 
                         minGSSize = 3, 
                         maxGSSize = 1000, 
                         pvalueCutoff = 0.05, 
                         verbose = TRUE, 
                         OrgDb = "org.Hs.eg.db", 
                         pAdjustMethod = "fdr",
                         seed = 1234)

#Make dataframe of gse results
gse_df <- gse@result %>% as.data.frame()
write.table(gse_df, "~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/GSEAiAMP21vsFusiondriven.tsv")
