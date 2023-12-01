## ---------------------------Title---------------------------------------------
# Comparing the model output to the shelf life data collected in the 2019 sampling 

## ---------------------------Description---------------------------------------
# Project: CIDA Spinach 

# Script description: See title

## ---------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(Metrics) 

## ---------------------------Data----------------------------------------------
#   Reading in raw data
model_output_inc.strains <- read.csv("inc.strains_total_count_by_package.csv", header = TRUE)
model_output_inc.temp <- read.csv("inc.temp_total_count_by_package.csv", header = TRUE)
model_output_inc.strains.temp <- read.csv("inc.strains.temp_total_count_by_package.csv", header = TRUE)
model_output_baseline <- read.csv("total_count_by_package.csv", header = TRUE)

## ---------------------------Checking model output distribution----------------
day <- 1:20
p_val <- vector(mode = "integer", length = length(day))

model_output_shapiro <- data.frame(cbind(day, p_val, p_val, p_val, p_val))

colnames(model_output_shapiro) <- c("day", "p_val_baseline", "p_val_inc.strains",
                            "p_val_inc_temp", "p_val_inc.strains.temp")

for(i in 1:20){
  base_data <- model_output_baseline$total_count[model_output_baseline$day == i]
  base_test <- shapiro.test(base_data)
  model_output_shapiro$p_val_baseline[model_output_shapiro$day == i] <- base_test$p.value
  
  inc.strains_data <- model_output_inc.strains$total_count[model_output_inc.strains$day == i]
  inc.strains_test <- shapiro.test(inc.strains_data)
  model_output_shapiro$p_val_inc.strains[model_output_shapiro$day == i] <- inc.strains_test$p.value
  
  inc.temp_data <- model_output_inc.temp$total_count[model_output_inc.temp$day == i]
  inc.temp_test <- shapiro.test(inc.temp_data)
  model_output_shapiro$p_val_inc_temp[model_output_shapiro$day == i] <- inc.temp_test$p.value
  
  inc.strains.temp_data <- model_output_inc.strains.temp$total_count[model_output_inc.strains.temp$day == i]
  inc.strains.temp_test <- shapiro.test(inc.strains.temp_data)
  model_output_shapiro$p_val_inc.strains.temp[model_output_shapiro$day == i] <- inc.strains.temp_test$p.value
}

#All data from the model output is not normally distributed. Will use the median observation for comparison with shelf life data 

med_shelf_life_baseline <- model_output_baseline %>%
  group_by(day) %>%
  summarize(median_baseline = median(total_count))

med_shelf_life_inc.strains <- model_output_inc.strains %>%
  group_by(day) %>%
  summarize(median_inc.strains = median(total_count))

med_shelf_life_inc.temp <- model_output_inc.temp %>%
  group_by(day) %>%
  summarize(median_inc.temp = median(total_count))

med_shelf_life_inc.strains.temp <- model_output_inc.strains.temp %>%
  group_by(day) %>%
  summarize(median_inc.strains.temp = median(total_count))

med_shelf_life_overall_p1 <- merge(med_shelf_life_baseline, med_shelf_life_inc.strains, by = "day")
med_shelf_life_overall_p2 <- merge(med_shelf_life_inc.temp, med_shelf_life_inc.strains.temp, by = "day")

med_shelf_life_overall <- merge(med_shelf_life_overall_p1, med_shelf_life_overall_p2)

## ------------------Ave. increase in microbial concentration-------------------

#Calculating the average increase in the median concentration of the packages, by: 
#(1) subtracting the microbial concentration of the packages in the baseline model from that of packages in the "trouble shooting" models,
#(3) averaging the median difference by the days used for this calculation

inc.strains_diff <- mean(med_shelf_life_overall$median_inc.strains - med_shelf_life_overall$median_baseline)
inc.strains_diff

inc.temp_diff <- mean(med_shelf_life_overall$median_inc.temp - med_shelf_life_overall$median_baseline)
inc.temp_diff

inc.strains.temp <- mean(med_shelf_life_overall$median_inc.strains.temp - med_shelf_life_overall$median_baseline)
inc.strains.temp

