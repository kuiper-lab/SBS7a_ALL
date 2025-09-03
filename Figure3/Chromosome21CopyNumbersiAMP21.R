#Load packages -----------------------------------------------------------------
library(data.table)
library(scales)
library(tidyverse)
library(cowplot)
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(readxl)
library(karyoploteR)

#Read in reference data --------------------------------------------------------
centro <- data.frame(fread('~/Projects/hypermutated_ALL/CODE/R/centromeres/cytoBand.hg38.centromeresOnly.txt'))
chromLength <- data.frame(fread('~/Downloads/chromosomeSizes.hg38.noCommas.txt'))
anonymous_names <- read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx")
metadata <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aPositiveiAMP21.xlsx")
metadata_biobank <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/SBS7aNegativeiAMP21.xlsx")

#Read in St Jude data ------------------------------------------------------------------
#SBS7a positive cases
SBS7apositive_iamp21_path_stjude <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/StJude/iamp21_SBS7aPositive/"
SBS7apositive_iamp21s_names_stjude <- c()
sample_idx <- 1
for (sample_path in list.files(SBS7apositive_iamp21_path_stjude, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), ".rdata")
  SBS7apositive_iamp21s_names_stjude <- c(SBS7apositive_iamp21s_names_stjude, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7apositive_iamp21s_stjude_gr <- sample_gr
  } else {
    SBS7apositive_iamp21s_stjude_gr <- c(SBS7apositive_iamp21s_stjude_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#SBS7a negative cases
SBS7anegative_iamp21_path_stjude <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/StJude/iamp21_SBS7aNegative/"
SBS7anegative_iamp21s_names_stjude <- c()

sample_idx <- 1
for (sample_path in list.files(SBS7anegative_iamp21_path_stjude, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), ".rdata")
  SBS7anegative_iamp21s_names_stjude <- c(SBS7anegative_iamp21s_names_stjude, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7anegative_iamp21s_stjude_gr <- sample_gr
  } else {
    SBS7anegative_iamp21s_stjude_gr <- c(SBS7anegative_iamp21s_stjude_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#Read in PMC data --------------------------------------------------------------
#SBS7a positive cases
SBS7apositive_iamp21_path_pmc <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/SBS7aPositiveiAMP21s/"
SBS7apositive_iamp21s_names_pmc <- c()
recode_list <- metadata$Anonymous_ID
names(recode_list) <- metadata$Tumor_ID
sample_idx <- 1
for (sample_path in list.files(SBS7apositive_iamp21_path_pmc, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), "_WGS.rdata")
  sample_name <- recode(sample_name,!!!recode_list)
  SBS7apositive_iamp21s_names_pmc <- c(SBS7apositive_iamp21s_names_pmc, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7apositive_iamp21s_pmc_gr <- sample_gr
  } else {
    SBS7apositive_iamp21s_pmc_gr <- c(SBS7apositive_iamp21s_pmc_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#SBS7a negative cases
SBS7anegative_iamp21_path_pmc <- "~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/granges/SBS7aNegativeiAMP21s/"
SBS7anegative_iamp21s_names_pmc <- c()
sample_idx <- 1
for (sample_path in list.files(SBS7anegative_iamp21_path_pmc, full.names = T)){
  load(sample_path)
  sample_name <- str_remove(basename(sample_path), "_.+")
  sample_name <- metadata_biobank$Anonymous_ID[grep(sample_name, metadata_biobank$Tumor_ID)]
  SBS7anegative_iamp21s_names_pmc <- c(SBS7anegative_iamp21s_names_pmc, sample_name)
  sample_gr$sample <- sample_name
  if (sample_idx == 1){
    SBS7anegative_iamp21s_pmc_gr <- sample_gr
  } else {
    SBS7anegative_iamp21s_pmc_gr <- c(SBS7anegative_iamp21s_pmc_gr, sample_gr)
  }
  sample_idx <- sample_idx + 1
}

#Plot coverage and include copy number -----------------------------------------
#Get coverage
chr21_weighted_coverage_sbs7apositive_stjude <- coverage(SBS7apositive_iamp21s_stjude_gr, weight = SBS7apositive_iamp21s_stjude_gr$COPY_NUMBER)[["chr21"]]
sbs7apositive_stjude_coverage <- setNames(c(runValue(chr21_weighted_coverage_sbs7apositive_stjude)/length(SBS7apositive_iamp21s_names_stjude),
                                            runValue(chr21_weighted_coverage_sbs7apositive_stjude)/length(SBS7apositive_iamp21s_names_stjude)),
                                          c(start(chr21_weighted_coverage_sbs7apositive_stjude),end(chr21_weighted_coverage_sbs7apositive_stjude)))
sbs7apositive_stjude_coverage <- sbs7apositive_stjude_coverage[as.character(sort(as.numeric(names(sbs7apositive_stjude_coverage))))]
chr21_weighted_coverage_sbs7anegative_stjude <- coverage(SBS7anegative_iamp21s_stjude_gr, weight = SBS7anegative_iamp21s_stjude_gr$COPY_NUMBER)[["chr21"]]
sbs7anegative_stjude_coverage <- setNames(c(runValue(chr21_weighted_coverage_sbs7anegative_stjude)/length(SBS7anegative_iamp21s_names_stjude),
                                            runValue(chr21_weighted_coverage_sbs7anegative_stjude)/length(SBS7anegative_iamp21s_names_stjude)),
                                          c(start(chr21_weighted_coverage_sbs7anegative_stjude),end(chr21_weighted_coverage_sbs7anegative_stjude)))
sbs7anegative_stjude_coverage <- sbs7anegative_stjude_coverage[as.character(sort(as.numeric(names(sbs7anegative_stjude_coverage))))]

chr21_weighted_coverage_sbs7apositive_pmc <- coverage(SBS7apositive_iamp21s_pmc_gr, weight = SBS7apositive_iamp21s_pmc_gr$COPY_NUMBER)[["chr21"]]
sbs7apositive_pmc_coverage <- setNames(c(runValue(chr21_weighted_coverage_sbs7apositive_pmc)/length(SBS7apositive_iamp21s_names_pmc),
                                            runValue(chr21_weighted_coverage_sbs7apositive_pmc)/length(SBS7apositive_iamp21s_names_pmc)),
                                          c(start(chr21_weighted_coverage_sbs7apositive_pmc),end(chr21_weighted_coverage_sbs7apositive_pmc)))
sbs7apositive_pmc_coverage <- sbs7apositive_pmc_coverage[as.character(sort(as.numeric(names(sbs7apositive_pmc_coverage))))]
chr21_weighted_coverage_sbs7anegative_pmc <- coverage(SBS7anegative_iamp21s_pmc_gr, weight = SBS7anegative_iamp21s_pmc_gr$COPY_NUMBER)[["chr21"]]
sbs7anegative_pmc_coverage <- setNames(c(runValue(chr21_weighted_coverage_sbs7anegative_pmc)/length(SBS7anegative_iamp21s_names_pmc),
                                         runValue(chr21_weighted_coverage_sbs7anegative_pmc)/length(SBS7anegative_iamp21s_names_pmc)),
                                       c(start(chr21_weighted_coverage_sbs7anegative_pmc),end(chr21_weighted_coverage_sbs7anegative_pmc)))
sbs7anegative_pmc_coverage <- sbs7anegative_pmc_coverage[as.character(sort(as.numeric(names(sbs7anegative_pmc_coverage))))]


#Set parameters
pp <- getDefaultPlotParams(plot.type = 1)
pp$data1min <- 0
pp$data1max <- 6

#Save plot
pdf("~/Projects/hypermutated_ALL/ANALYSES/cnvAnalysis/iamp21_mean_copynumber.pdf",
    width = 10)

#Make base
kp <- plotKaryotype(chromosomes=c("chr21"), plot.params = pp, genome = "hg38")

#Plot
kpLines(kp,
        chr = "chr21",
        x = as.numeric(names(sbs7apositive_stjude_coverage)),
        y = sbs7apositive_stjude_coverage,
        col = "orange")
kpLines(kp,
        chr = "chr21",
        x = as.numeric(names(sbs7anegative_stjude_coverage)),
        y = sbs7anegative_stjude_coverage,
        col = "black")
kpLines(kp,
        chr = "chr21",
        x = as.numeric(names(sbs7apositive_pmc_coverage)),
        y = sbs7apositive_pmc_coverage,
        col = "yellow")
kpLines(kp,
        chr = "chr21",
        x = as.numeric(names(sbs7anegative_pmc_coverage)),
        y = sbs7anegative_pmc_coverage,
        col = "lightgrey")
kpRect(kp,
       chr = "chr21",
       y0 = 0, y1 = 6,
       x0 = 20e6, x1 = 25e6,
       border = "black", col = NA)
kpRect(kp,
       chr = "chr21",
       y0 = 0, y1 = 6,
       x0 = 40e6, x1 = 46709983,
       border = "black", col = NA)
kpAxis(kp)
dev.off()
