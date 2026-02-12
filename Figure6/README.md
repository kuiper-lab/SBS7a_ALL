# Figure 6
This directory contains code to analyse single-cell whole genome sequencing data and generate the lineage tree for patient P0625. These scripts are also used for the generation of supplementary figures S9 and S10.

## VAFMQFilter.R
This script can be used to apply stricter filters to the filtered vcfs resulting from PTATO v1.3.3. In this case, we select mutations with at least 5 supporting alternative reads, an allele frequency of at least 0.25 and a mapping quality of at least 59.

## GetHaplotypeCallerPositions.sh
This script uses the unfiltered vcf from PTATO v1.3.3 and the sample-specific filtered vcf files generated with VAFMQFilter.R to generate one filtered vcf file that contains all positions which will be considered for variant calling.

## generateHaplotypeCallerJob.sh
This script creates jobs to run haplotypecaller using the positions generated with GetHaplotypeCallerPositions.sh

## generateGenotypeJob.sh
This script combines the gvcfs generated with generateHaplotypeCallerJob.sh and genotypes them. This results in one vcf with calls for the single cells, the MSCs and the tumor bulk sample. This can be used as input for CellPhyWrapper.

## P0625LineageTree.R
This script uses the output of CellPhyWrapper to generate a lineage tree. It extracts mutational signatures from each branch and plots the contributions on the tree (Figure 6A). It also plots the mutational profiles and allele frequencies for SBS7a-positive branches (Figure 6B&C). Additionally, it plots some branch-specific quality metrics (Supplementary Figure S9).

## fit_to_signatures_strict_tree_method.R and plot_tree_contributions_branchnames.R
These scripts contain functions which are used in P0625LineageTree.R.
