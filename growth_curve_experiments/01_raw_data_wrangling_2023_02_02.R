  ## -------------------------------------------------------------------------
# Wrangling Raw Data for Growth Curve Experiments Conducted at 6C

## -------------------------------------------------------------------------
# Project: CIDA Spinach 

# Script description: combining the raw data to a single data frame. Ensuring that all variables are in of the appropriate type. Wrangling data as necessary (e.g., log transformation). 

## -------------------------------------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(lubridate)

## -------------------------------------------------------------------------
#   Reading in raw data
b2 <- read.csv("data/raw/batch2_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b3 <- read.csv("data/raw/batch3_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b4 <- read.csv("data/raw/batch4_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b5 <- read.csv("data/raw/batch5_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b6 <- read.csv("data/raw/batch6_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b7 <- read.csv("data/raw/batch7_consolidated_data_2022_1_8.csv", header = TRUE, stringsAsFactors = TRUE)
b8 <- read.csv("data/raw/batch8_consolidated_data_2022_05_17.csv", header = TRUE, stringsAsFactors = TRUE)
b9 <- read.csv("data/raw/batch9_consolidated_data_2022_05_17.csv", header = TRUE, stringsAsFactors = TRUE)

## -------------------------------------------------------------------------
# Merging data into a single dataframe; Formatting variables;
combined_data <- dplyr::bind_rows(b2, b3, b4, b5, b6, b7, b8, b9) #All dataframes need have the same column names to use this function. NAs are added when column names do not between dataframes

#Check whether any NAs were introduced during merging
cbind(
  lapply(
    lapply(combined_data, is.na)
    , sum))

#Select variables necessary for analysis
combined_data <- combined_data %>%
  select(c(Plate.Id, Counted.Colonies, Concentration, Dilution.Factor))

#Making the column names snake case
colnames(combined_data) <- gsub("[.]", "_", colnames(combined_data)) %>%
  tolower() 

#Separating the Plate.Id column, to provide more resolution between each observation. Plate.Id is the name entered for each observation in SphereFlash. It can be broken down into descriptive variables (e.g., day, media etc)
combined_data <- combined_data %>%
separate(plate_id, into = c("batch", "temperature", "day", "isolate", "biorep", "media", "dilution", "time", "platingrep"), sep = "_", remove = FALSE) 

#Check whether any NAs were introduced during separating Plate.ID
cbind(
  lapply(
    lapply(combined_data, is.na)
    , sum))

#Convert combined_data$Batch into a factor variable
combined_data$batch <- as.factor(combined_data$batch)

#Convert combined_data$Day, combined_data$BioRep, combined_data$Dilution, combined_data$PlatingRep, combined_data$Concentration into a factor variable
combined_data$day <- as.numeric(combined_data$day)
combined_data$biorep <- as.numeric(combined_data$biorep)
combined_data$dilution <- as.numeric(combined_data$dilution)
combined_data$platingrep <- as.numeric(combined_data$platingrep)
combined_data$concentration <- as.numeric(combined_data$concentration)

#Renaming the concentration vector, to indicate that this concentration was calculated by SphereFlash. 
colnames(combined_data)[12] <- "concentration_from_sphereflash"

#End of merging data and formatting variables 

## -------------------------------------------------------------------------
# Separating data frames for further analysis; Additional formatting;

#Separate the negative controls and the spinach controls, as a separate data set
control_data <- combined_data %>%
  filter(isolate == "NC" | isolate == "Control")

#Check whether any NAs were introduced during filtering out negative control and spinach control data
cbind(
  lapply(
    lapply(control_data, is.na)
    , sum))

#Separate the APC and strain growth data 
apc_strain_data <- combined_data %>%
  filter(isolate != "NC" & isolate != "Control")

#Check whether any NAs were introduced during filtering out APC and strain growth data 
cbind(
  lapply(
    lapply(apc_strain_data, is.na)
    , sum))

#End of separating data and additional formatting 

## -------------------------------------------------------------------------
# Save wrangled data to R Project folder 

#Push the wrangled data back to the R project 
date <- Sys.Date()
date <- gsub("-", "_", date)

write.csv(apc_strain_data, paste("data/wrangled/", date, "_apc_strain_data_wrangled_01.csv", sep = ""), row.names = FALSE)
write.csv(control_data, paste("data/wrangled/", date, "_control_data_wrangled_01.csv", sep = ""), row.names = FALSE)

# End of saving data into R Project folder 

## -------------------------------------------------------------------------
# Data frame Guide:

#   b2 - b9: Data collected from each individual experiment. Each batch consists of three independent/biological replicates. 
#   combined_data <- Data from b2 - b9 merged into a single data frame




