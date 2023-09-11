## -----------------------------Title-------------------------------------------
# Wrangling Raw Data for Growth Curve Experiments Conducted at 6C

## -----------------------------Description-------------------------------------
# Project: CIDA Spinach 

# Script description: Wrangling data after parsing in OpenRefine. This wrangling is conducted to address plates that had 0 CFU during shelf life testing and also to re-structure the data for future analysis.

## -----------------------------Packages----------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); 

## -----------------------------Reading in Data---------------------------------
apc_strain_data <- read.csv("data/wrangled/2023_02_02_apc_strain_data_wrangled_01.csv", header = TRUE)

## -----------------------------Wrangling---------------------------------------

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
#25% of the detection limit
apc_strain_data_wrangled_ave$log_average_wrangled_conc <-
  ifelse(apc_strain_data_wrangled_ave$average_wrangled_conc == 0,
         log10(0.25*apc_strain_data_wrangled_ave$dilution_factor/0.05),
         apc_strain_data_wrangled_ave$log_average_wrangled_conc)
  
#Removing unecessary columns from the apc_strain_data dataframe
apc_strain_data_wrangled_ave <- apc_strain_data_wrangled_ave %>%
  select(-c("concentration_from_sphereflash", "counted_colonies", "dilution", "dilution_factor", "platingrep"))
 
#Identify which day the max count is reached for each isolate. Counts will be marked as death_phase or not 
max_by_isolate_and_rep <- apc_strain_data_wrangled_ave %>%
  group_by(batch, isolate, media, time) %>%
  arrange(desc(log_average_wrangled_conc)) %>%
  slice(1)

## -------------------------------------------------------------------------
#Save wrangled data to R Project folder 

#Push the wrangled data back to the R project 
date <- Sys.Date()
date <- gsub("-", "_", date)

write.csv(apc_strain_data_wrangled_ave, paste("data/wrangled/", date, "_apc_strain_averaged_wrangled_02.csv", sep = ""), row.names = FALSE)
write.csv(max_by_isolate_and_rep, paste("data/wrangled/", date, "_max_counts_for_dieoff_assessment_wrangled_02.csv", sep = ""), row.names = FALSE)
# End of saving data into R Project folder 

