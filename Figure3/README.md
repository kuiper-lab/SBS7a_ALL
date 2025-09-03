# Figure 3
This directory contains code to plot karyotypes, compare prognosis based on karyotype, and study copy number alterations of chromosome 21 in more detail. CNAs were called according to GATK best practices, using GATK v4.1.7.0 and an in-house panel of normals.

## PlotHyperdiploidKaryotypes.R
This script takes manually assessed karyotypes based on copy number profiles and plots and clusters them (Figure 3A). It also plots the percentage of tumors with a full chromosome gain for each chromosome (Figure 3A). Additionally, it stratifies each sample into a good or poor prognosis based on the karyotype, and plots this per SBS7a status (Figure 3B).

## Chromosome21CopyNumbersiAMP21.R
This script can be used to calculate the average copy number across chromosome 21, using segment files generated with GATK (Figure 3D).

## RegionsChromosome21iAMP21.R
This script can be used to zoom in to certain regions of chromosome 21 and plot the average copy number based on segment files generated with GATK. The script also calculates the average copy number across each region per sample and plots these values (Figure 3D).
