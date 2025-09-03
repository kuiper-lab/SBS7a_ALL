# Figure S6
This directory contains code to analyse single-cell whole genome sequencing data and generate the lineage tree for patient P0624. It uses the same steps as for patient P0625, for which the code can be found in directory Figure6.

## P0624LineageTree.R
This script uses the output of CellPhyWrapper to generate a lineage tree. It extracts mutational signatures from each branch and plots the contributions on the tree (Figure S6A). It also plots the mutational profiles and allele frequencies for SBS7a-positive branches (Figure S6B,S6C). Additionally, it plots some branch-specific quality metrics (Figure S7).