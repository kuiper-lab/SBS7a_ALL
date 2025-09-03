# Figure 4
This directory contains code to identify subclonal mutations from somatic vcfs generated with Mutect2 from GATK v4.1.1.0 and extract mutational signatures from these mutations. These scripts are also used for the generation of supplementary figures S4 and S5.

## AnalyseSubclonalMutations.R
This script is used to identify samples with mutations of which the allele frequency follows a bimodal distribution, and clusters the mutations into a clonal and subclonal fraction by using the antimode as a threshold. In this process, violin plots of the allele frequency are created (Supplementary Figure S4). For tumors that show a bimodal distribution with at least 200 subclonal mutations, mutational signatures are extracted for the clonal and subclonal mutations (Figure 4A&B, Supplementary Figure S5A&B). The input for this script are somatic vcfs generated with Mutect2 from GATK v4.1.1.0.

## AnalyseSubclonalMutationsWithOverlap.R
This script can be used to apply stricter threshold, when the allele frequencies of the clonal and subclonal mutations show great overlap. It uses manual threshold to identify the clonal and subclonal mutations, and extracts mutational signatures from both groups (Figure 4C, Supplementary Figure S5C). The input for this script are somatic vcfs generated with Mutect2 from GATK v4.1.1.0.

## CombineFiguresSubclonalMutations.R
This script is used to combine the sample specific plots into one figure (Figure 4, Supplementary Figure S5). It uses the plots generated with AnalyseSubclonalMutations.R and AnalyseSubclonalMutationsWithOverlap.R.