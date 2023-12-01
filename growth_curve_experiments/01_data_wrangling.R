## -----------------------------Title-------------------------------------------
# Wrangling Raw Data for Growth Curve Experiments Conducted at 6C, after marking counts in death phase on OpenRefine

## -----------------------------Description-------------------------------------
# Project: CIDA Spinach 

# Script description: Wrangling data after parsing in OpenRefine. 
#This wrangling is conducted to remove data points that were in the death phase. 

#Additionally, calculating the SM-APC of the spinach

## -----------------------------Packages----------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); 

## -----------------------------Reading in Data---------------------------------
apc_strain_data <- read.csv("data/apc_strain_averaged_wrangled_01_openrefine_deathphase.csv", header = TRUE)

## -----------------------------Wrangling---------------------------------------
#Subtracting strain-specific counts from the BHI counts
#Using data from 36h timepoint for S12-0166, S12-0180 and S12-0184
#Using data from 24h for S12-0116, S12-0132, S12-0144
apc_strain_data_2 <- apc_strain_data %>%
  filter((isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "36") |
           (isolate %in% c("S12-0116", "S12-0132", "S12-0141"))) 

#There was not BHI observation for Batch 7 (B7), S12-0184, day 2
#Thus, SM-APC (in the apc_strain_data_bhi_adj) will not be calculated for this data point
apc_strain_data_bhi_adj <- apc_strain_data_2 %>%
  filter(!(isolate == "S12-0184" & batch == "B7" & day == "2" & temperature == "6")) %>%
  group_by(temperature, batch, day, isolate) %>%
  mutate(average_wrangled_conc_adj = average_wrangled_conc[media == "BHI"] - average_wrangled_conc[media == "BHIrif"])

apc_strain_data_filtered_p1 <- apc_strain_data_bhi_adj %>%
  filter(media == "BHI") %>%
  filter(death_phase == "no") %>%
  select(-(average_wrangled_conc))

colnames(apc_strain_data_filtered_p1)[11] <- "average_wrangled_conc"
  
#Removing counts that are in the death phase. Keep counts that from 36h counts only for S12-0166, S12-0180 and S12-0184
apc_strain_data_filtered_p2 <- apc_strain_data %>%
  filter(death_phase == "no") %>%
  filter((isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "36") |
  (isolate %in% c("S12-0116", "S12-0132", "S12-0141"))) %>%
  filter(media == "BHIrif")

apc_strain_filtered <- bind_rows(apc_strain_data_filtered_p1, apc_strain_data_filtered_p2)

#Removing SM-APC data, which was equal to or below 0  CFU/g
apc_strain_final <- apc_strain_filtered %>%
  filter(!(media == "BHI" & average_wrangled_conc <= 0))

#Calculating the log-transformed concentration for the SM-APC observations

apc_strain_final$log_average_wrangled_conc <- ifelse(apc_strain_final$media == "BHI", 
                                                     log10(apc_strain_final$average_wrangled_conc),
                                                     apc_strain_final$log_average_wrangled_conc)

apc_strain_death_phase <- apc_strain_data %>%
  filter(death_phase == "yes") 

apc_strain_24h_slow_growing_isolates <- apc_strain_data %>%
  filter(death_phase == "no") %>%
  filter(isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "24")

apc_level_lower_than_strain <- apc_strain_filtered %>%
  filter(average_wrangled_conc < 0)

## -----------------------------Export------------------------------------------
#Save wrangled data to R Project folder 

#Push the wrangled data back to the R project 
write.csv(apc_strain_final, paste("data/apc_strain_averaged_wrangled_02.csv", sep = ""), row.names = FALSE)

