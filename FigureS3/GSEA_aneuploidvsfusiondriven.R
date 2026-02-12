#General---------------------------------------------------------------------
#Load libraries
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(readxl)
library(cowplot)

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
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "fdr",
             seed = 1234)

#Make dataframe of GO terms
gse_full_df <- gse@result %>% as.data.frame()

#Get list of GO terms that were significant in both separate analyses
intersect_df <- read.table("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/hyperdiploid_intersect.tsv")

#Plot
pdf("~/Projects/hypermutated_ALL/ANALYSES/RNAseq/DEanalysis_09022026/gse_intersect_dotplot.pdf", width = 12)
dotplot(gse[gse@result$ID %in% intersect_df$ID, asis = TRUE], 
        x = "GeneRatio",
        showCategory = 50, # Choose how many categories should be showed
        font.size = 12, 
        split = ".sign",
        title = "Enriched GO terms Aneuploid vs Fusion-driven") + 
  facet_grid(.sign ~ ., scales = "free", space = "free") +
  scale_y_discrete(labels = function(gse) str_wrap(gse, width = 100)) +
  theme(strip.text.y = element_text(size = 13),
        plot.title = element_text(face = "bold", size = 14))
dev.off()
