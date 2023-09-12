## -------------------------------------------------------------------------
# Wrangling Raw Data for Growth Curve Experiment on Baby Spinach 
## -------------------------------------------------------------------------
# Project: CIDA Spinach 

# Script description: combining the raw data to a single data frame. 
#Ensuring that all variables are in of the appropriate type. 
#Wrangling data as necessary (e.g., log transformation). 

## -------------------------------------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(lubridate)

## -------------------------------------------------------------------------
#   Reading in raw data
b2 <- read.csv("data/raw/batch2_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b3 <- read.csv("data/raw/batch3_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b4 <- read.csv("data/raw/batch4_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b5 <- read.csv("data/raw/batch5_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b6 <- read.csv("data/raw/batch6_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b7 <- read.csv("data/raw/batch7_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b8 <- read.csv("data/raw/batch8_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)
b9 <- read.csv("data/raw/batch9_consolidated_data.csv", header = TRUE, stringsAsFactors = TRUE)

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

apc_strain_data[, c("batch", "temperature", "isolate", "biorep", "media", "time")] <- 
  lapply(apc_strain_data[, c("batch", "temperature", "isolate", "biorep", "media", "time")] , as.factor)

#Calculating the average concentration for each day of testing and media. Removing columns that would not be required during subsequent portions of the analysis
apc_strain_data_wrangled_ave <- apc_strain_data %>%
  group_by(batch, isolate, day, time, media) %>%
  mutate(average_wrangled_conc = mean(concentration_from_sphereflash)) %>%
  distinct(batch, isolate, day, time, media, .keep_all = TRUE) %>%
  mutate(log_average_wrangled_conc= log10(average_wrangled_conc)) %>%
  ungroup()

#Identifying which observations are in below the detection limit 
below_det_limit <- apc_strain_data_wrangled_ave[apc_strain_data_wrangled_ave$average_wrangled_conc == 0, ]

#Replacing the concentration of samples that were below the detection limit with 
#25% of the detection limit (detection limit = dilution factor*volume plated (50uL))
apc_strain_data_wrangled_ave$log_average_wrangled_conc <-
  ifelse(apc_strain_data_wrangled_ave$average_wrangled_conc == 0,
         log10(0.25*apc_strain_data_wrangled_ave$dilution_factor/0.05),
         apc_strain_data_wrangled_ave$log_average_wrangled_conc)

#Removing unecessary columns from the apc_strain_data dataframe
apc_strain_data_wrangled_ave <- apc_strain_data_wrangled_ave %>%
  select(-c("concentration_from_sphereflash", "counted_colonies", "dilution", "dilution_factor", "platingrep"))

#Identify which day the max count is reached for each isolate. Counts will be marked as death_phase or not, based on this data 
max_by_isolate_and_rep <- apc_strain_data_wrangled_ave %>%
  group_by(batch, isolate, media, time) %>%
  arrange(desc(log_average_wrangled_conc)) %>%
  slice(1)

## -------------------------------------------------------------------------
#Save wrangled data to R Project folder 

#Push the wrangled data back to the R project 

write.csv(apc_strain_data_wrangled_ave, "data/wrangled/apc_strain_averaged_wrangled_01.csv", row.names = FALSE)
write.csv(max_by_isolate_and_rep, "data/wrangled/max_counts_for_dieoff_assessment_wrangled_01.csv", row.names = FALSE)
# End of saving data into R Project folder 

