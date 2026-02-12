# Supplementary Figure S5
This directory contains code to study the timing of SBS7a compared to amplifications of chromosome 21. CNAs were called according to GATK best practices, using GATK v4.1.7.0 and an in-house panel of normals.

## FilterCNACalls.R
This script can be used to extract only CNAs from segment files generated with GATK.

## FilterVcfsForCNAAnalysis.R
This script filters somatic SBS vcfs generated with Mutect2 from GATK v4.1.1.0 in a similar way as for the other analyses, but keeps mutations with an allele frequency of at least 0.1.

## GetSomaticGainedSBSs.R
This script selects somatic mutations from the filtered vcf which overlap with an amplification of chromosome 21, using the output from FilterCNACalls.R and FilterVcfsForCNAAnalysis.R.

## GetGermlineGainedSBSs.R
This script selects germline variants which are heterozygous in the germline (allele frequency between 0.25 and 0.75) and present in the tumor (allele frequency of at least 0.1). From these variants it selects the ones which overlap with amplifications of chromosome 21. It used the output from FilterCNACalls.R and GATK HaplotypeCaller.

## MakeSNVCNATimingPlot.R
This script can be used to plot the allele frequency of germline and somatic mutations to study the timing of SBS7a compared to amplifications of chromosome 21 (Supplementary Figure S5). It used the output from GetSomaticGainedSBSs.R and GetGermlineGainedSBSs.R.
