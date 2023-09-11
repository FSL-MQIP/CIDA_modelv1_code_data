## ----------------------Title--------------------------------------------------
#   Summarizing and plotting the plating data collected during the 2019 sampling

## ----------------------Description--------------------------------------------
#   Project: CIDA Spinach 

#  Script description: Summarizing the change in APC, GN PC, Y and M over processing and shelf life for the 2019 sampling 

## ----------------------Packages-----------------------------------------------
library(tidyverse); library(ggplot2); library(dplyr)

## ------------------------------Raw Data---------------------------------------
#Reading in the raw data 

spinachmass <- read.csv("data/raw/0919A_Sample_Mass.csv", header = TRUE) 
inprocess <- read.csv("data/raw/PlateCountData_PilotSampling_InProcessSamples_Sept2019_101019_parsed_2022_08_23.csv", header = TRUE)
shelflife <- read.csv("data/raw/PlateCountData_PilotSampling_ShelfLifeSamples_Sept2019_102619_parsed_2022_08_23.csv", header = TRUE)

## ------------------------------Wrangling--------------------------------------

#Subsetting sample mass data by in-process or shelf life
samplemass_inprocess <- spinachmass[1:22, ]
samplemass_shelflife <- spinachmass[23:55, ]

#For sample mass data, splitting the sampleID vector into unique elements/vectors (sampling, sample_id, rep)
samplemass_inprocess <- samplemass_inprocess %>%
  separate(col = sampleID, into = c("sampling", "sample_id", "rep"), sep = "_")

samplemass_shelflife <- samplemass_shelflife %>%
  separate(col = sampleID, into = c("sampling", "rep", "sample_id"), sep = "_")

#If the spinach was weighed directly into the bag, calculating the total mass of the bag by adding the sample mass with the buffer mass
#The total bag mass will be used for downstream analysis 
samplemass_inprocess$total_mass <- ifelse(is.na(samplemass_inprocess$total_mass), 
                                          samplemass_inprocess$sample_mass + samplemass_inprocess$buffer_mass, 
                                          samplemass_inprocess$total_mass)

samplemass_shelflife$total_mass <- ifelse(is.na(samplemass_shelflife$total_mass), 
                                          samplemass_shelflife$sample_mass + samplemass_shelflife$buffer_mass, 
                                          samplemass_shelflife$total_mass)

#For the microbial count data, separated the sampleID vector into unique elements/vectors (sampling, sample_id, rep)
#Created a "df" column, to store the dilution factor of the enumerated plates 
inprocess <- inprocess %>%
  separate(col = sampleID, into = c("sampling", "sample_id", "rep"), sep = "-") %>%
  select(-c("date_plated", "date_enumerated", "photoID", "notes")) %>%
  filter(!is.na(raw_count)) %>%
  mutate(df = NA)

shelflife <- shelflife %>%
  separate(col = sampleID, into = c("sampling", "rep", "sample_id"), sep = "_") %>%
  select(-c("date_plated", "date_enumerated", "photoID", "notes", "X")) %>%
  filter(!grepl("BHI", test)) %>%
  filter(!is.na(raw_count)) %>%
  mutate(df = NA)

#Checking whether there are any observations where the microbial count was below the detection limit 

#Only one plate replicate was below the detection limit: sample_id "D", replicate "B", test "M", plate rep "1". 
#No need for substituting the detection limit, since the other plate replicate was above the detection limit (as data from
#plate replicates will be averaged). 
table(inprocess$raw_count == "0")
table(shelflife$raw_count == "0")

#Checking whether both plate replicate were enumerated from the same dilution
#If so, the raw counts can be averaged and then the microbial concentration can be calculated
inprocess %>%
  group_by(sample_id, rep, test) %>%
  mutate(same_dilution = dilution[1] == dilution[2]) %>%
  with(table(same_dilution))

shelflife %>%
  group_by(sample_id, rep, test) %>%
  mutate(same_dilution = dilution[1] == dilution[2]) %>%
  with(table(same_dilution))

#Identifying which of the enumerated plate replicates from the shelf life data set 
#were not from the same dilution in the shelf life data set
shelflife_check <- shelflife %>%
  group_by(sample_id, rep, test) %>%
  mutate(same_dilution = dilution[1] == dilution[2]) 

shelflife[(shelflife_check$same_dilution == "FALSE"), ]

#Since there are some plate replicates that were from different dilutions, the 
#microbial concentration will be calculated for each plate replicate and then averaged 

#Defining a function for calculating the dilution factor 
dilution_factor_func <- function(x, y, z){
  #x = sample mass, y = sample mass + buffer mass (i.e., total mass), z = dilution of enumerated plate
  first_dilution = log10(y/x)
  dilution_plated = z
  dilution_factor = 10^(first_dilution + (dilution_plated - 1))
  return(dilution_factor)
}

#Calculating the dilution factor for each plate replicate 
for(i in 1:nrow(inprocess)){
  row_index <- which(inprocess$sample_id[i] == samplemass_inprocess$sample_id & inprocess$rep[i] == samplemass_inprocess$rep)
  inprocess$df[i] <- dilution_factor_func(samplemass_inprocess$sample_mass[row_index], samplemass_inprocess$total_mass[row_index], inprocess$dilution[i])
}

for(i in 1:nrow(shelflife)){
  row_index <- which(shelflife$sample_id[i] == samplemass_shelflife$sample_id & shelflife$rep[i] == samplemass_shelflife$rep)
  shelflife$df[i] <- dilution_factor_func(samplemass_shelflife$sample_mass[row_index], samplemass_shelflife$total_mass[row_index], shelflife$dilution[i])
}

#Calculating the microbial concentration for each sample, by multiplying the raw count with the dilution factor
inprocess <- inprocess %>%
  mutate(conc = raw_count * df)

shelflife <- shelflife %>%
  mutate(conc = raw_count * df)

#Averaging the microbial concentration for each sample replicate, by microbiological test. I.e., the microbial concentrations 
#that were calculated based on plate replicates are averaged, to obtain a single concentration for each sample replicate 
#and microbiological test
inprocess_final <- inprocess %>%
  group_by(sample_id, rep, test) %>%
  mutate(average_conc = mean(conc)) %>%
  distinct(sample_id, rep, test, .keep_all = TRUE) %>%
  select(-c("plate_rep", "df", "dilution", "raw_count", "conc")) %>%
  filter(!grepl("NS", rep))

shelflife_final <- shelflife %>%
  group_by(sample_id, rep, test) %>%
  mutate(average_conc = mean(conc)) %>%
  distinct(sample_id, rep, test, .keep_all = TRUE) %>%
  select(-c("plate_rep", "df", "dilution", "raw_count", "conc"))

## ------------------------------Export-----------------------------------------

wrangled_counts <- bind_rows(inprocess_final, shelflife_final)

date <- Sys.Date()
date <- gsub("-", "_", date)

#write.csv(wrangled_counts, paste0("data/wrangled/0919a_microbial_counts_", date, ".csv"), row.names = FALSE)
