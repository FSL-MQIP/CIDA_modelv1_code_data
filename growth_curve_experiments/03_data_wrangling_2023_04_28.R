## -----------------------------Title-------------------------------------------
# Wrangling Raw Data for Growth Curve Experiments Conducted at 6C, after marking counts in death phase on OpenRefine

## -----------------------------Description-------------------------------------
# Project: CIDA Spinach 

# Script description: Wrangling data after parsing in OpenRefine. This wrangling is conducted to address plates that had 0 CFU during shelf life testing and also to re-structure the data for future analysis.

## -----------------------------Packages----------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); 

## -----------------------------Reading in Data---------------------------------
apc_strain_data <- read.csv("data/wrangled/2023_02_03_apc_strain_averaged_wrangled_02_openrefine_deathphase.csv", header = TRUE)

## -----------------------------Wrangling---------------------------------------
#Subtracting strain-specific counts from the BHI counts
apc_strain_data_2 <- apc_strain_data %>%
  filter((isolate %in% c("S12-0166", "S12-0180", "S12-0184") & time == "36") |
           (isolate %in% c("S12-0116", "S12-0132", "S12-0141"))) 

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

apc_strain_final <- apc_strain_filtered %>%
  filter(average_wrangled_conc > 0)

apc_strain_final$log_average_wrangled_conc <- log10(apc_strain_final$average_wrangled_conc)

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
date <- Sys.Date()
date <- gsub("-", "_", date)

write.csv(apc_strain_filtered, paste("data/wrangled/", date, "_apc_strain_averaged_wrangled_03.csv", sep = ""), row.names = FALSE)
write.csv(apc_strain_death_phase, paste("data/wrangled/", date, "_apc_strain_counts_in_death_phase_wrangled_03.csv", sep = ""), row.names = FALSE)
write.csv(apc_strain_24h_slow_growing_isolates, paste("data/wrangled/", date, "_apc_strain_24h_counts_slow_growing_isolates_wrangled_03.csv", sep = ""), row.names = FALSE)
write.csv(apc_level_lower_than_strain, paste("data/wrangled/", date, "_apc_level_lower_than_strain_wrangled_03.csv", sep = ""), row.names = FALSE)

