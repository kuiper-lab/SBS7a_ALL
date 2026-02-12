# Figure S3
This directory contains code to analyse differences in transcription between BCP-ALL subtypes, in this case comparing aneuploid BCP-ALL (iamp21 and high hyperdiploid) to ETV6::RUNX1 BCP-ALL.

## iAMP21vsETV6RUNX1.R and HyperdiploidvsETV6RUNX1.R
These scripts use RNAseq data from 58 cases of ETV6::RUNX1 BCP-ALL and either 10 cases of iAMP21 of 72 cases of high hyperdiploid BCP-ALL to find differentially expressed genes and gene sets.

## GetIntersectiAMP21Hyperdiploid.R
This script finds overlapping genes and gene sets which are differentially expressed in iAMP21 and high hyperdiploid BCP-ALL compared to ETV6::RUNX1 BCP-ALL. It performs over-representation analysis on the differentially expressed genes and plots these (Figure S4).

## GSEA_aneuploidvsfusiondriven.R
This script performs differential expression analysis and gene set enrichment analysis comparing iAMP21 and high hyperdiploid BCP-ALL to ETV6::RUNX1 BCP-ALL, to make plots of the overlapping gene sets identified with GetIntersectiAMP21Hyperdiploid.R.
