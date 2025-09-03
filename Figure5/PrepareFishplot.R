#Load packages
library(tidyverse)

#Load cluster vaf data
P0608_vafs <- read.table("~/Projects/hypermutated_ALL/ANALYSES/SBS7a/clustered_SNVs_nogermlinefilter_nosoftclippedbases_23042025//P0608-clustered_SNVs.txt",
                         row.names = 1)

#Rename clusters
P0608_vafs$cluster[P0608_vafs$cluster == "4"] <- "ClR1"
P0608_vafs$cluster[P0608_vafs$cluster == "12b"] <- "ClDx"
P0608_vafs$cluster[P0608_vafs$cluster == "2"] <- "ClR2"
P0608_vafs$cluster[P0608_vafs$cluster == "5"] <- "ClDx_fall"

#Get mean vaf per cluster
for (cluster in unique(P0608_vafs$cluster)){
  mean_vaf_Dx <- mean(P0608_vafs$XP0608Dx_read_based_AF[P0608_vafs$cluster == cluster])
  mean_vaf_R1 <- mean(P0608_vafs$XP0608R1_read_based_AF[P0608_vafs$cluster == cluster])
  mean_vaf_R2 <- mean(P0608_vafs$XP0608R2_read_based_AF[P0608_vafs$cluster == cluster])
  mean_vaf_R3 <- mean(P0608_vafs$XP0608R3_read_based_AF[P0608_vafs$cluster == cluster])
  if (cluster == unique(P0608_vafs$cluster)[1]){
    cluster_df <- data.frame(cluster = cluster, Dx = mean_vaf_Dx, R1 = mean_vaf_R1, R2 = mean_vaf_R2, R3 = mean_vaf_R3)
  } else {
    cluster_df <- rbind(cluster_df, c(cluster, mean_vaf_Dx, mean_vaf_R1, mean_vaf_R2, mean_vaf_R3))
  }
}

#Calculate CCF
cluster_df$CCF_Dx <- as.numeric(cluster_df$Dx) * 2
cluster_df$CCF_R1 <- as.numeric(cluster_df$R1) * 2
cluster_df$CCF_R2 <- as.numeric(cluster_df$R2) * 2
cluster_df$CCF_R3 <- as.numeric(cluster_df$R3) * 2

#Save
write.csv(cluster_df, "~/Projects/hypermutated_ALL/ANALYSES/deepsequencing/CCF_P0608_06052025.csv")
