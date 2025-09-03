# Figure 5
This directory contains code to study acquired mutations, analyse deep sequencing data and generate fishplots.

## ClusterSBS7aMutations.R
This script uses somatic vcfs generated with Mutect2 from GATK v4.1.1.0, and filters and clusters SBSs based on the dynamics of their allele frequencies over time.

## RefitClusteredMutations.R
This script is used to extract mutational signatures from the clusters of mutations generated with ClusterSBS7aMutations.R.

## AcquiredMutationsPlot.R
This script can be used to plot the acquired mutations for each relapse based on the contributions calculated with RefitClusteredMutations.R (Figure 5).

# RunBamReadCount.sh
This script can be used to analyze the raw alternative read count for deep sequencing data with bam-readcount, using reads mapped to the GRCh38 human reference genome with BWA v0.7.13.

# DeepSequencing.R
This script was used to analyse the deep sequencing data based on the files generated with RunBamReadCount.sh. It also plots the allele frequency for each mutation (Figure 5).

## PrepareFishplot.R
This script can be used to transform the format of clustered mutations generated with ClusterSBS7aMutations.R so they can be used to make a fishplot.

## Fishplot.R
This script can be used to make fishplots using the output from PrepareFishplot.R (Figure 5).