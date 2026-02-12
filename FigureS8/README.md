# Figure S8
This directory contains code to find pathogenic driver genes and calculate the probability that they were caused by SBS7a.

## MakepathogenicSNVlist.R
This script uses clustered vcfs to create a table of mutations, which are filtered on CADD score (at least 15) and a list of known driver genes.

## GetSNVprobabilities.R
This script takes the mutation list from MakepathogenicSNVlist.R, and uses signature contributions to calculate a probability of being caused by SBS7a for each mutation.

## PlotPathogenicSNVs.R
This script does some additional filtering and plots the remaining mutations.
