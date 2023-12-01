## ----------------------Title--------------------------------------------------
#   Summarizing and plotting the plating data collected during the 2019 sampling

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

#  Script description: Summarizing and plotting the change in APC, GN PC, Y and M over processing and shelf life for the 2019 sampling 

## ----------------------Packages-----------------------------------------------
library(tidyverse); library(dplyr); library(stringr)

## --------------------------------Data-----------------------------------------
#Reading in the data 
data <- read.csv("data/0919a_microbial_counts.csv", header = TRUE)

##-------------------------------Wrangling--------------------------------------
#Taking the geometric average and mean of the sample replicates 
data_log <- data %>%
  mutate(log_conc = log10(average_conc)) %>%
  dplyr::select(-average_conc)

#Separating the inprocess and shelf life data 

inprocess_samples <- c("H", "T", "V", "S", "W", "D")

inprocess <- data_log %>%
  filter(sample_id %in% inprocess_samples)

shelflife <- data_log %>%
  filter(!(sample_id %in% inprocess_samples))


#Summarizing the average microbial count at different stages of the supply chain

summary_inprocess <- inprocess %>%
  group_by(test, sample_id) %>%
  summarize(average = round(mean(log_conc), 2), sd = round(sd(log_conc), 2))

high_inprocess <- summary_inprocess %>%
  group_by(test) %>%
  arrange(desc(average)) %>%
  slice(1)

low_inprocess <- summary_inprocess %>%
  group_by(test) %>%
  arrange(average) %>%
  slice(1)

#Reduction after washing

wash_reduction <- summary_inprocess %>%
  group_by(test) %>%
  summarize(wash_reduction = average[sample_id == "S"] - average[sample_id == "W"])

#Identifying the highest microbial count observed during shelf life, by test 

summary_shelflife <- shelflife %>%
  group_by(test, sample_id) %>%
  summarize(average = round(mean(log_conc), 2), sd = round(sd(log_conc), 2)) %>%
  mutate(day = substr(sample_id, start = 2, stop = nchar(sample_id)))


high_shelflife <- summary_shelflife %>%
  group_by(test) %>%
  arrange(desc(average)) %>%
  slice(1)

#Calculating the absolute increase in microbial counts by test 
summary_abs_shelflife <- summary_shelflife %>%
  group_by(test) %>%
  summarize(absolute_increase = max(average) - min(average))


##-------------------------------Export-----------------------------------------

write.csv(summary_inprocess, "outputs/summary_inprocess.csv")
write.csv(summary_shelflife, "outputs/summary_shelflife.csv")


  
