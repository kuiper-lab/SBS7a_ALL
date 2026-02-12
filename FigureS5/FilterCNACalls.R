#Load packages
library(tidyverse)

#Get called segment file
segment_files <- list.files("~/Projects/hypermutated_ALL/DATA/cnv/UV_20082024/",
                            full.names = T, pattern = ".called.seg", recursive = T)

#Get tumor files
segment_files <- c(grep("tumor", segment_files, value = T))

#Loop over files
for (segment_file in segment_files){
  #Get data from file
  con = file(segment_file, "r")
  segment_data = readLines(con)
  close(con)
  sample <- str_remove(basename(segment_file), "\\..+\\.called\\.seg")
  
  #Get length of header
  header_length = length(grep("@|#", substr(segment_data,1,1)))
  
  #Read in data
  segment_df <- read_tsv(segment_file, skip = header_length)
  
  #Take only called regions
  CNAs <- segment_df[segment_df$CALL != 0,]
  
  #Save file with CNAs
  write_tsv(CNAs, paste0("~/Projects/hypermutated_ALL/RESULTS/cnvAnnotated/GATK_calls/", sample, ".tsv"))
}
