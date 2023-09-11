##--------------------------Title-----------------------------------------------
# 0919a isolates: Analyzing blast results and assigning metadata

##--------------------------Description-----------------------------------------
#The matches in blast with the highest percent identity will be selected for each
#isolate. The top match will be manually checked, to observed that the coverage is 
#not low (> 95% coverage of the query sequence). Each isolate will be linked with
#it's previous ID, to connect the isolates with information about the sample of origin
#and the day of isolation. 

##--------------------------Packages--------------------------------------------
library(dplyr); library(tidyverse)

##--------------------------Raw Data--------------------------------------------
blast <- read.csv("data/raw/0919a_raw_blast_results.csv", header = TRUE)
metadata <- read.csv("data/raw/FMT_isolates.csv", header = TRUE)

##--------------------------Wrangling-------------------------------------------
top_match <- blast %>%
  group_by(qseqid) %>%
  arrange(desc(pident)) %>%
  slice(1)

top_match_pruned <- top_match %>%
  select(qseqid, subject.titles) 

top_match_pruned$qseqid <- gsub("0919a", "", top_match_pruned$qseqid)

metadata_pruned <- metadata %>%
  select(FSL, Previous.ID)
  
fsl_previd <- merge(top_match_pruned, metadata_pruned, by.x = "qseqid", by.y = ("FSL"), all = FALSE)


##--------------------------Export----------------------------------------------

date <- Sys.Date()
date <- gsub("-", "_", date)

#write.csv(fsl_previd, paste0("data/wrangled/", date, "_0919a_isolates_fsl_previd.csv"), row.names = FALSE)
