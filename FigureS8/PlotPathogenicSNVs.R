#Load packages
library(ggplot2)
library(tidyverse)

#Load color palette
Colorpalette<-scale_fill_manual(values = c("SBS1" = "#9BC0CD",
                                           "SBS5" = "#A9A9A9",
                                           "SBS7a" = "#F6EB13",
                                           "SBS2" = "#CE2627",
                                           "SBS13" = "#B12325",
                                           "SBS18" = "#AA6AAC",
                                           "SBSA" = "#DADADA",
                                           "SBS86" = "#0A0F25",
                                           "SBS87" = "#5C6BC0",
                                           "SBS17a" = "#FFC2D1",
                                           "SBS17b" = "#FF8FAB"))

#Read file with pathogenic variants
Pathogenic_SNVs<-read.delim(file="20251217_SNVprobabilities_pathogenic.tsv",header=T)

# Add gene annotation
for (i in 1:nrow(Pathogenic_SNVs)){
  Pathogenic_SNVs$Patient_gene_amino_acid_change[i]<-paste(Pathogenic_SNVs$anonymous[i],
                                                           Pathogenic_SNVs$gene[i],
                                                           Pathogenic_SNVs$protein_position[i],
                                                           Pathogenic_SNVs$amino_acid[i])
  
}

# Select patients with over 50% SBS7a contribution in a mutation, 24 mutations left
Pathogenic_SNVs_SBS7a <- Pathogenic_SNVs[Pathogenic_SNVs$SBS7a >= 0.50,]

# Select for match with SBS7a trinucleotide context, 16 mutations left
Pathogenic_SNVs_SBS7a <- Pathogenic_SNVs_SBS7a[Pathogenic_SNVs_SBS7a$Trinucleotide_context_matches_SBS7a == "Yes",]

# Filter on CADD 20 instead of CADD15. This removes numerous non exonic variants. 12 variants left, of which 4 intronic.
Pathogenic_SNVs_SBS7a <- Pathogenic_SNVs_SBS7a[Pathogenic_SNVs_SBS7a$CADD >= 20,]

# Prepare file for plotting
Pathogenic_SNVs_SBS7a_long<-gather(Pathogenic_SNVs_SBS7a,Signature,Probability,SBSA:SBS87)
Pathogenic_SNVs_SBS7a_long$Probability[Pathogenic_SNVs_SBS7a_long$Probability == 0]<-NA

#Plot SBS7a variants barplot
Probability_barplot_SBS7a_final<- ggplot(data=Pathogenic_SNVs_SBS7a_long,
                                         aes(y=factor(Patient_gene_amino_acid_change, levels = c("P0813 KRAS 12 G/S",
                                                                                                 "P0611 ZEB2 NA NA",
                                                                                                 "P0625 CREBBP 768 R/*",
                                                                                                 "P0610 ELL NA NA",
                                                                                                 "P0631 CREBBP NA NA",
                                                                                                 "P0612 KRAS 146 A/T",
                                                                                                 "P0608 NR3C1 612 S/L",
                                                                                                 "P0624 CREBBP 1490 D/N",
                                                                                                 "P0145 KMT2D 203 S/F",
                                                                                                 "P0608 NR3C1 386 R/*",
                                                                                                 "P0617 IKZF1 507 S/L",
                                                                                                 "P0608 KDM6A NA NA")),
                                             x=Probability,
                                             fill=factor(Signature,levels = c("SBS87","SBS7a")))) +
  geom_bar(position="stack", stat = "identity", width = 0.8) +
  Colorpalette +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour="black",linewidth =1.5),
        axis.text.x = element_text(angle= 60, hjust = 1),
        axis.ticks.y = element_line(colour = "black",linewidth=1.5)) +
  labs(x=NULL,fill="Mutational\nSignature", y = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02), labels = scales::percent)
