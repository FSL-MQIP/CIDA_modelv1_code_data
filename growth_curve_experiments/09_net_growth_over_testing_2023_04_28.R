##  ----------------------------Title-------------------------------------------
#Assessing the net change in inoculated strains and APC over shelf life
##  --------------------------Description---------------------------------------
#   Project: CIDA Project 

#   Script description: See title
##  --------------------------Packages------------------------------------------
library(dplyr)
library(reshape2)

##----------------------------Data----------------------------------------------
apc_strain_data <- read.csv("data/wrangled/2023_04_28_apc_strain_averaged_wrangled_03.csv", header = TRUE)

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
  summarize(max_conc = mean(growth_mean), max_conc_sd = sd(growth_mean))

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
  summarize(max_conc = mean(growth_mean), max_conc_sd = sd(growth_mean))

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

apc_conc_max_6 <- apc_strain_data %>%
  filter(media == "BHI" & temperature == 6) %>%
  group_by(temperature, isolate, biorep, day) %>%
  summarize(growth_mean = round(mean((log_average_wrangled_conc)), 2), 
            growth_sd = round(sd((log_average_wrangled_conc)), 2)) %>%
  group_by(isolate, biorep) %>% 
  arrange(desc(growth_mean)) %>%
  slice(1) %>% 
  ungroup() %>%
  group_by(biorep) %>%
  summarize(max_conc = mean(growth_mean), max_conc_sd = sd(growth_mean),
            min_conc = mean(growth_mean), min_conc_sd = sd(growth_mean))

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

