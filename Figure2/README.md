# Figure 2
This directory contains code for the calculation and visualization of tumor mutational load, number of CC>TT mutations and transcriptional strand bias.

## GetSBS7aPositiveSkinCancer.R
This script is used to extract mutational signatures from filtered vcf files of adult skin cancers, received from the Hartwig Medical Foundation, to select the SBS7a-positive cases.

## MakeMutationalLoadPlot.R
This script is used to plot the tumor mutational load for SBS7a-negative BCP-ALL, SBS7a-positive BCP-ALL, SBS7a-positive ALCL and SBS7a-positive skin cancer (Figure 2A). This is done based on the generated SBS mutation matrices.

## GetMNPmatrix.R
This script is used to filter somatic vcfs containing multiple nucleotide polymorphisms, to generate DBS mutational matrices. The somatic vcfs were either generated with Mutect2 from GATK v4.1.1.0 or received from the Hartwig Medical Foundation in the case of adult skin cancers.

## MakeDBS1Plot.R
This script is used to plot the number of CC>TT mutations compared to the number of C>T mutations (Figure 2B). It takes the mutation matrices generated with GetMNPmatrix.R as input.

## CalculateTranscriptionalStrandBias.R
This script filters vcf files generated with Mutect2 from GATK v4.1.1.0 and calculates transcriptional strand bias.

## PlotTSB.R
This script is used to plot the transcriptional strand bias for SBS7a-negative BCP-ALL, SBS7a-positive BCP-ALL, SBS7a-positive ALCL and SBS7a-positive skin cancer (Figure 2C). The input for this script is generated with CalculateTranscriptionalStrandBias.R.
