#Load packages
library(tidyverse)
library(readxl)

#Load data
absolute_contributions_SBS7a_pertimepoint <- read.delim("~/Projects/hypermutated_ALL/RESULTS/tables/AbsoluteContributions_SBS7aALLNoClusters_SBS7aPerTimepoint.tsv", sep = "\t")
metadata <- read_excel("~/Projects/hypermutated_ALL/DOCS/SBS7aCohort_metadata.xlsx")
anonymous_names <- read_excel("~/surfdrive/Shared/Kuiper group/Relapsed_ALL/DOCS/ALL_Relapse_Coded_Num.xlsx")

#Sort on subtype
mutational_burdens <- colSums(absolute_contributions_SBS7a_pertimepoint)
mutational_burdens_unique <- c()
for (patient in str_extract(names(mutational_burdens), "P\\d\\d\\d\\d+")){
  highest_burden <- names(which.max(mutational_burdens[grep(patient, names(mutational_burdens))]))
  mutational_burdens_unique[patient] <- mutational_burdens[highest_burden]
}
cohort_overview <- as.data.frame(mutational_burdens_unique)
cohort_overview$patient <- rownames(cohort_overview)
colnames(cohort_overview) <- c("burden", "patient")
cohort_overview$subtype <- NA
for (patient in cohort_overview$patient){
  print(patient)
  subtype <- metadata$`Genetic subtype`[grep(patient, metadata$`Anonymous ID`)]
  cohort_overview[patient, "subtype"] <- subtype
}
cohort_overview$subtype[cohort_overview$subtype == "Down Syndrome/Hyperdiploid"] <- "High Hyperdiploid"
cohort_overview$subtype <- factor(cohort_overview$subtype,
                                  levels = c("High Hyperdiploid",
                                             "iAMP21",
                                             "Low Hypodiploid",
                                             "B-Other"))
cohort_overview_sorted <- arrange(cohort_overview, subtype, burden)
patient_order <- rownames(cohort_overview_sorted)
for (patient in patient_order){
  patient_df <- as.data.frame(absolute_contributions_SBS7a_pertimepoint[, sort(grep(patient, colnames(absolute_contributions_SBS7a_pertimepoint), value = T))])
  colnames(patient_df) <- sort(grep(patient, colnames(absolute_contributions_SBS7a_pertimepoint), value = T))
  if (patient == patient_order[1]){
    absolute_contributions_SBS7a_pertimepoint_sorted <- patient_df
  } else {
    absolute_contributions_SBS7a_pertimepoint_sorted <- cbind(absolute_contributions_SBS7a_pertimepoint_sorted, patient_df)
  }
}

#Transform table
absolute_contributions_rotated <- as.data.frame(t(absolute_contributions_SBS7a_pertimepoint_sorted))
absolute_contributions_rotated$Sample <- rownames(absolute_contributions_rotated)
absolute_contributions_long <- pivot_longer(absolute_contributions_rotated,
                                            grep("SBS", colnames(absolute_contributions_rotated), value = T),
                                            names_to = "Signatures",
                                            values_to = "SigLoad")

#Add anonymous names
absolute_contributions_long$Patient <- str_extract(absolute_contributions_long$Sample, "P\\d\\d\\d\\d+")

#Add white lines
patients <- unique(str_extract(absolute_contributions_long$Sample, "P\\d\\d\\d\\d+"))
sample_levels <- c()
for(patient_idx in 1:(length(patients)-1)){
  patient <- patients[patient_idx]
  last_row <- max(grep(patient, absolute_contributions_long$Sample))
  empty_sample <- data.frame(Sample = paste0("blank", patient_idx), Signatures = unique(absolute_contributions_long$Signatures), SigLoad = 0, Patient = NA)
  absolute_contributions_long <- rbind(absolute_contributions_long[1:last_row, ], empty_sample, absolute_contributions_long[(last_row+1):nrow(absolute_contributions_long), ])
  sample_levels <- c(sample_levels, unique(grep(patient, absolute_contributions_long$Sample, value = T)), paste0("blank", patient_idx))
}
sample_levels <- c(sample_levels, unique(grep(patients[length(patients)], absolute_contributions_long$Sample, value = T)))

absolute_contributions_long$Sample <- factor(absolute_contributions_long$Sample, levels = sample_levels)

#Put SBS7a at the bottom
absolute_contributions_long$Signatures <- factor(absolute_contributions_long$Signatures,
                                                 levels = c(unique(absolute_contributions_long$Signatures[absolute_contributions_long$Signatures != "SBS7a"]), "SBS7a"))

#Plot
Colorpalette <- scale_fill_manual(values = c("SBSA"= "#DADADA",
                                             "SBS1" = "#9BC0CD",
                                             "SBS2" = "#CE2627",
                                             "SBS7a" = "#F6EB13",
                                             "SBS13" = "#B12325",
                                             "SBS18" = "#AA6AAC",
                                             "SBS87" = "#5C6BC0"))


overviewplot<-ggplot(data=absolute_contributions_long,aes(x=Sample,y=SigLoad,fill=Signatures)) +
  geom_col(position= "stack",colour="white") +
  labs(y="Number of SNVs",x=NULL) + 
  Colorpalette +
  scale_x_discrete(name = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.ticks.x = element_blank(), 
        aspect.ratio = 3/5,
        panel.background = element_blank(),
        legend.position = c(0.1,0.7), 
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm')) 
overviewplot

pdf("~/Projects/hypermutated_ALL/RESULTS/plots/paper/CohortOverview.pdf",
    width = 15,
    height = 10)
overviewplot +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        axis.ticks.x = element_blank(), 
        aspect.ratio = 3/5,
        panel.background = element_blank(),
        legend.position = c(0.05,0.8), 
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm')) 
dev.off()

#Make subtype plot
Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")
absolute_contributions_long$Subtype <- NA
absolute_contributions_subtypes <- absolute_contributions_long[grep("blank", absolute_contributions_long$Sample, invert = T), ]
for (patient in unique(str_extract(absolute_contributions_subtypes$Sample, "P\\d\\d\\d\\d+"))){
  print(patient)
  subtype <- as.character(cohort_overview$subtype[grep(patient, cohort_overview$patient)])
  absolute_contributions_subtypes$Subtype[grep(patient, absolute_contributions_subtypes$Sample)] <- subtype
}

absolute_contributions_long[grep("blank", absolute_contributions_long$Sample, invert = T), "Subtype"] <- absolute_contributions_subtypes$Subtype
absolute_contributions_long[grep("blank", absolute_contributions_long$Sample), "Subtype"] <- "B-Other"

subtype_plot <- ggplot(data=absolute_contributions_long,aes(x=Sample,y=SigLoad,fill=Subtype)) +
  geom_col(position= "stack",colour="white") +
  labs(y="Number of SNVs",x=NULL) + 
  scale_fill_manual(values = c("lightgrey", Tol_light[c(1,3,5,4)])) +
  scale_x_discrete(name = element_blank(),
                   labels=element_blank()) +
  theme(axis.text.x = element_text(face="bold"),
        axis.ticks.x = element_blank(), 
        aspect.ratio = 3/5,
        panel.background = element_blank(),
        legend.position = c(0.2,0.8), 
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm'))
subtype_plot

pdf("~/Projects/hypermutated_ALL/RESULTS/plots/paper/CohortOverview_subtypes.pdf",
    width = 15,
    height = 10)
subtype_plot
dev.off()

