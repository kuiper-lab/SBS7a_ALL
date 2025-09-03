#load packages
library(fishplot)
library(tidyverse)

#Load data
cluster_df <- read.csv("~/Projects/hypermutated_ALL/ANALYSES/deepsequencing/CCF_P0608_06052025.csv")

timepoints <- c(0,1,2,3,4,5,6,7,8,9)
parents <- c(0,1,1,3)
fracTable <- matrix(c(cluster_df[c(2,3,1,4), c("CCF_Dx")],
                      0.01,0,0,0,
                      0.01,0,0,0,
                      cluster_df[c(2,3,1,4), c("CCF_R1")],
                      0.01,0,0.009,0,
                      0.01,0,0.009,0,
                      cluster_df[c(2,3,1,4), c("CCF_R2")],
                      0.01,0,0.009,0.008,
                      0.01,0,0.009,0.008,
                      cluster_df[c(2,3,1,4), c("CCF_R3")]),
                    ncol = length(timepoints))
fracTable <- fracTable * 100
fracTable[2,1] <- 89

#Create fish object
fish = createFishObject(fracTable,
                        parents = parents,
                        timepoints = timepoints)

#calculate the layout of the drawing
fish = layoutClones(fish)

#draw the plot, using the splining method (recommended)
#and providing both timepoints to label and a plot title
pdf("~/Projects/hypermutated_ALL/ANALYSES/fishplots/deepseq/P0608_07052025.pdf", width = 10, height = 5)
fishPlot(fish,
         shape="spline",
         title.btm="P0608",
         vlines=c(0,3,6,9), 
         vlab=c("Initial Diagnosis","Relapse 1", "Relapse 2", "Relapse 3"),
         bg.type = "solid",
         bg.col = "white",
         col.vline = "black",
         cex.title = 2,
         cex.vlab = 2)
dev.off()
