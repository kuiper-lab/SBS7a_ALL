# Figure 1
This directory contains code for the generation of mutation matrices and extraction of mutational signatures. It also contains code for the visualizations in figure 1.

## FilterSBS7aPatients.R
This script is used to filter somatic vcfs generated with Mutect2 from GATK v4.1.1.0. It generates filtered vcfs and mutation matrices.

## DeNovoExtraction.R
This script is used to extract mutational signatures from the mutation matrices generated with FilterSBS7aPatients.R.

## RefitSignatures.R
This script is used to select signatures based on the de novo extraction with DeNovoExtraction.R, and perform a refit using the selected signatures. It outputs both relative and absolute contributions.

## CohortOverview.R
This script is used to plot an overview of the absolute contributions of the cohort as generated with RefitSignatures.R (Figure 1A).

## SBS7aPercentagePerSubtypeAndAge.R
This script is used to plot the percentage of SBS7a-positive tumors per subtype as determined with RefitSignatures.R (Figure 1B). It also plots the age at diagnosis for SBS7a-positive and SBS7a-negative tumors per subtype (Figure 1C).