#Load packages -----------------------------------------------------------------
source("createdFunctions_HM-ALL.R") #From https://github.com/kuiper-lab/MultipleRelapse
library(MutationalPatterns)
library(tidyverse)
library(cowplot)
library(NMF)
library(readxl)
library(ggsignif)
cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#Load metadata
metadata_subtype <- read_excel("~/Projects/hypermutated_ALL/ANALYSES/metadata/B-ALL_subtypes_manual.xlsx")
metadata_SBS7a <- read_excel("~/Projects/hypermutated_ALL/DOCS/SBS7aCohort_metadata.xlsx")

#Load mutation matrix
load("~/Projects/hypermutated_ALL/ANALYSES/deNovoExtraction/rdata/SBS7aNoClusters_MutationMatrix.rdata")

#Make dataframe with subtypes and SBS7a-status
overview_df <- metadata_subtype
overview_df <- overview_df[overview_df$Anonymous %in% str_extract(colnames(sbs_mutation_matrix_combined), "P\\d+"),]
overview_df$SBS7aStatus <- "SBS7a-negative"
overview_df$SBS7aStatus[overview_df$Anonymous %in% metadata_SBS7a$`Anonymous ID`] <- "SBS7a-positive"

#Get SBS7a status per subtype
subtype_counts <- as.data.frame(table(overview_df[, c("FinalSubtype", "SBS7aStatus")]))
colnames(subtype_counts) <- c("Subtype", "SBS7a Status", "Number of Samples")
subtype_counts$Subtype <- factor(subtype_counts$Subtype, levels = names(sort(table(overview_df$FinalSubtype), decreasing = T)))

#Contract small subtypes to "Other Subtype"
subtype_numbers <- table(overview_df$FinalSubtype)
other_subtypes <- names(subtype_numbers)[subtype_numbers == 1]
other_subtypes_positive <- sum(subtype_counts$`Number of Samples`[subtype_counts$Subtype %in% other_subtypes &
                                                                    subtype_counts$`SBS7a Status` == "SBS7a-positive"])
other_subtypes_negative <- sum(subtype_counts$`Number of Samples`[subtype_counts$Subtype %in% other_subtypes &
                                                                    subtype_counts$`SBS7a Status` == "SBS7a-negative"])
subtype_counts_short <- subtype_counts[!subtype_counts$Subtype %in% other_subtypes,]
subtype_counts_short <- rbind(subtype_counts_short,
                                       data.frame("Subtype" = "Other Subtype",
                                                  "SBS7a Status" = c("SBS7a-positive",
                                                                     "SBS7a-negative"),
                                                  "Number of Samples" = c(other_subtypes_positive,
                                                                          other_subtypes_negative),
                                                  check.names = F))

#Simplify and add numbers
subtype_counts_short_simplified <- subtype_counts_short
subtype_counts_short_simplified$n <- table(overview_df$FinalSubtype)[as.character(subtype_counts_short_simplified$Subtype)]
subtype_counts_short_simplified$n[(nrow(subtype_counts_short_simplified)-1):nrow(subtype_counts_short_simplified)] <-
  other_subtypes_positive + other_subtypes_negative
subtype_counts_short_simplified$Subtype <- str_remove(subtype_counts_short_simplified$Subtype, " B-ALL")
subtype_counts_short_simplified$Subtype <- paste0(subtype_counts_short_simplified$Subtype,
                                                           " (n=",
                                                           subtype_counts_short_simplified$n,
                                                           ")")
subtype_counts_short_simplified_percentage <-
  data.frame(
    subtypes = unique(subtype_counts_short_simplified$Subtype),
    percentage = subtype_counts_short_simplified$`Number of Samples`[subtype_counts_short_simplified$`SBS7a Status` == "SBS7a-positive"] /
      (
        subtype_counts_short_simplified$`Number of Samples`[subtype_counts_short_simplified$`SBS7a Status` == "SBS7a-positive"] +
          subtype_counts_short_simplified$`Number of Samples`[subtype_counts_short_simplified$`SBS7a Status` == "SBS7a-negative"]
      )
  )
recode_list <- subtype_counts_short_simplified$n
names(recode_list) <- subtype_counts_short_simplified$Subtype
subtype_counts_short_simplified_percentage$n <- recode(subtype_counts_short_simplified_percentage$subtypes, !!!recode_list)
subtype_counts_short_simplified$Subtype <- factor(subtype_counts_short_simplified$Subtype,
                                                           levels = arrange(subtype_counts_short_simplified_percentage, desc(percentage), desc(n))$subtypes)

#Plot
subtype_plot_simplified <- ggplot(subtype_counts_short_simplified) +
  geom_col(aes(x = Subtype, y = `Number of Samples`, fill = `SBS7a Status`), position = "fill") +
  scale_fill_manual(values = c("lightgrey", cbp2[4])) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.position = "right",
        strip.text = element_text(size = 15)) +
  ggtitle("Occurrence Of SBS7a In An\nUnselected BCP-ALL Cohort (n=191)") +
  ylab("Fraction Of Patients") +
  xlab("BCP-ALL Subtype")
ggsave(plot = subtype_plot_simplified,
       filename = "~/Projects/hypermutated_ALL/RESULTS/plots/paper/SBS7aPercentagePerSubtype.pdf",
       width = 7.5, height = 6, units = "in")

#Plot age only for subtypes where SBS7a is seen
overview_df_selectedsubtypes <- overview_df[overview_df$FinalSubtype %in% c("High Hyperdiploid B-ALL",
                                                                            "iAMP21 B-ALL",
                                                                            "Low Hypodiploid B-ALL"),]
overview_df_selectedsubtypes$Subtype <- str_remove(overview_df_selectedsubtypes$FinalSubtype,
                                                   " B-ALL")
ggplot(overview_df_selectedsubtypes) +
  geom_boxplot(aes(x = Subtype,
                   y = `Age at diagnosis`,
                   fill = SBS7aStatus)) +
  theme_classic() +
  scale_fill_manual(values=c("lightgrey", cbp2[4])) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        title = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        legend.position = "right",
        strip.text = element_text(size = 15)) +
  ggtitle("Age At Diagnosis Of Patients With BCP-ALL") +
  ylab("Age At Diagnosis") +
  xlab("BCP-ALL Subtype")
ggsave(filename = "~/Projects/hypermutated_ALL/RESULTS/plots/paper/SBS7avsAge.pdf",
       width = 7.5, height = 6, units = "in")

model <- lm(`Age at diagnosis` ~ SBS7aStatus + FinalSubtype + SBS7aStatus * FinalSubtype,
             data = overview_df_selectedsubtypes)
drop1(model, test = "F")
model2 <- lm(`Age at diagnosis` ~ SBS7aStatus + FinalSubtype,
             data = overview_df_selectedsubtypes)
drop1(model2, test = "F")