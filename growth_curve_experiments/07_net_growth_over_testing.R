##  ----------------------------Title-------------------------------------------
#Assessing the net change in inoculated strains and APC over shelf life
##  --------------------------Description---------------------------------------
#   Project: CIDA Project 

#   Script description: See title

#For the SM-APC data, the lots are labelled as follows:
#Lot 1, 2, 3 correspond to biorep 2, 3, 4 of FSL S12-0116, FSL S12-0132 and FSL S12-0141, at 6˚C
#Lot 4, 5, 6 correspond to biorep 1, 2, 3 of FSL S12-0166, FSL S12-0180 and FSL S12-0184, at 6˚C
#Lot 7 corresponds to biorep 1 of FSL S12-0116, FSL S12-0132 and FSL S12-0141, at 10˚C
#Lot 7 corresponds to biorep 1 of FSL S12-0166, FSL S12-0180 and FSL S12-0184, at 10˚C

##  --------------------------Packages------------------------------------------
library(dplyr)
library(reshape2)

##----------------------------Data----------------------------------------------
apc_strain_data <- read.csv("data/apc_strain_averaged_wrangled_02.csv", header = TRUE)

##----------------------------Data Wrangling------------------------------------

mutant_strain_net_growth <- apc_strain_data %>%
  filter(media == "BHIrif") %>%
  group_by(temperature, isolate, biorep) %>%
  summarize(net_growth = max(log_average_wrangled_conc) - min(log_average_wrangled_conc)) %>%
  group_by(temperature, isolate) %>%
  summarize(mean_net_growth = round(mean(net_growth), 2), sd_net_growth = round(sd(net_growth), 2))

ave_mutant_strain_net_growth <- mutant_strain_net_growth %>%
  group_by(temperature) %>%
  summarize(mean = round(mean(mean_net_growth), 2))

mutant_conc_d0 <- apc_strain_data %>%
  filter(media == "BHIrif" & day == 0) %>%
  group_by(temperature, isolate) %>%
  summarize(d0_growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            d0_growth_sd = round(sd((log_average_wrangled_conc)), 2))

mutant_conc_max_6 <- apc_strain_data %>%
  filter(media == "BHIrif" & temperature == 6) %>%
  group_by(temperature, isolate, biorep, day) %>%
  summarize(growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            growth_sd = round(sd((log_average_wrangled_conc)), 2)) %>%
  group_by(isolate, biorep) %>% 
  arrange(desc(growth_mean)) %>%
  slice(1) %>% 
  ungroup() %>%
  group_by(isolate) %>%
  summarize(max_conc = round(mean(growth_mean), 2), max_conc_sd = round(sd(growth_mean), 2))

mutant_conc_max_10 <- apc_strain_data %>%
  filter(media == "BHIrif" & temperature == 10) %>%
  group_by(temperature, isolate, biorep, day) %>%
  summarize(growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            growth_sd = round(sd((log_average_wrangled_conc)), 2)) %>%
  group_by(isolate, biorep) %>% 
  arrange(desc(growth_mean)) %>%
  slice(1) %>% 
  ungroup() %>%
  group_by(isolate) %>%
  summarize(max_conc = round(mean(growth_mean),2), max_conc_sd = round(sd(growth_mean), 2))

apc_net_growth_1 <- apc_strain_data %>%
  filter(media == "BHI") %>%
  filter(isolate == "S12-0116" | isolate == "S12-0132" | isolate == "S12-0141") %>%
  group_by(temperature, isolate, biorep) %>%
  summarize(net_growth = max(log_average_wrangled_conc) - min(log_average_wrangled_conc)) %>%
  group_by(temperature, biorep) %>%
  summarize(mean_net_growth = round(mean(net_growth), 2), sd_net_growth = round(sd(net_growth), 2))

apc_net_growth_2 <- apc_strain_data %>%
  filter(media == "BHI") %>%
  filter(isolate == "S12-0166" | isolate == "S12-0180" | isolate == "S12-0184") %>%
  group_by(temperature, isolate, biorep) %>%
  summarize(net_growth = max(log_average_wrangled_conc) - min(log_average_wrangled_conc)) %>%
  group_by(temperature, biorep) %>%
  summarize(mean_net_growth = round(mean(net_growth), 2), sd_net_growth = round(sd(net_growth), 2))

apc_net_growth <- rbind(apc_net_growth_1, apc_net_growth_2)

ave_apc_net_growth <- apc_net_growth %>%
  group_by(temperature) %>%
  summarize(mean = round(mean(mean_net_growth), 2))

apc_conc_d0 <- apc_strain_data %>%
  filter(media == "BHI" & day == 0) %>%
  group_by(temperature, isolate) %>%
  summarize(d0_growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            d0_growth_sd = round(sd((log_average_wrangled_conc)), 2))

apc_max_conc_1 <- apc_strain_data %>%
  filter(media == "BHI") %>%
  filter(isolate == "S12-0116" | isolate == "S12-0132" | isolate == "S12-0141") %>%
  group_by(temperature, isolate, biorep) %>%
  summarize(max_biorep = max(log_average_wrangled_conc)) %>%
  group_by(temperature, biorep) %>%
  summarize(max = round(mean(max_biorep), 2))

apc_max_conc_2 <- apc_strain_data %>%
  filter(media == "BHI") %>%
  filter(isolate == "S12-0166" | isolate == "S12-0180" | isolate == "S12-0184") %>%
  group_by(temperature, isolate, biorep) %>%
  summarize(max_biorep = max(log_average_wrangled_conc)) %>%
  group_by(temperature, biorep) %>%
  summarize(max = round(mean(max_biorep), 2))

apc_max_conc_6 <- rbind(apc_max_conc_1, apc_max_conc_2)

apc_conc_max_10 <- apc_strain_data %>%
  filter(media == "BHI" & temperature == 10) %>%
  group_by(temperature, isolate, biorep, day) %>%
  summarize(growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            growth_sd = round(sd((log_average_wrangled_conc)), 2)) %>%
  group_by(isolate, biorep) %>% 
  arrange(desc(growth_mean)) %>%
  slice(1) %>% 
  ungroup() %>%
  group_by(isolate) %>%
  summarize(max_conc = mean(growth_mean), max_conc_sd = sd(growth_mean))

