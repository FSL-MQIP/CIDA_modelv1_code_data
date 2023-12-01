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
model_output <- read.csv("total_count_by_package.csv", header = TRUE)
sampling_data <- read.csv("0919a_microbial_counts.csv", header = TRUE)
packages_summary <- read.csv("packages_summary.csv", header = TRUE)
  
## ---------------------------Checking model output distribution----------------
day <- 1:20
p_val <- vector(mode = "integer", length = length(day))

model_output_shapiro <- data.frame(cbind(day, p_val))

for(i in 1:20){
  data <- model_output$total_count[model_output$day == i]
  test <- shapiro.test(data)
  model_output_shapiro$p_val[model_output_shapiro$day == i] <- test$p.value
}

#All data from the model output is not normally distributed. Will use the median observation for comaprison with shelf life data 

med_shelf_life <- model_output %>%
  group_by(day) %>%
  summarize(median = median(total_count))

## ---------------------------Data Wrangling------------------------------------

shelf_life_data <- sampling_data %>%
  filter(grepl("D", sample_id)) %>%
  filter(!(sample_id %in% c("D", "D0")) & test == "APC")

shelf_life_data$sample_id <- gsub("D", "", shelf_life_data$sample_id)

shelf_life_data$sample_id <- as.numeric(shelf_life_data$sample_id)

med_shelf_life <- med_shelf_life %>%
  filter(day %in% shelf_life_data$sample_id)

## ---------------------------RSME----------------------------------------------

#Calculating the geometrically averaged mean observation for each day of sampling 
shelf_life_data_rmse <-shelf_life_data %>%
  mutate(log_conc = log10(average_conc)) %>%
  group_by(sampling, sample_id) %>%
  summarize(geom_log_conc = mean(log_conc))

rmse_auto <- rmse(shelf_life_data_rmse$geom_log_conc, med_shelf_life$median)
rmse_auto

## ---------------------------Accuracy factor-----------------------------------
#Calculating the geometrically averaged mean observation for each day of sampling 
shelf_life_data_af <-shelf_life_data_rmse %>%
  mutate(geom_ln_conc = log(geom_log_conc)) %>%
  dplyr::select(-geom_log_conc)

med_shelf_life_af <- med_shelf_life %>%
  mutate(median_ln = log(median)) %>%
  dplyr::select(-median)

af_sum <- 0

for(i in 1:nrow(shelf_life_data_af)){
  d <- shelf_life_data_af$sample_id[i]
  add <- (med_shelf_life_af$median_ln[med_shelf_life_af$day == d] -
            shelf_life_data_af$geom_ln_conc[shelf_life_data_af$sample_id == d])^2
  af_sum <- af_sum + add
}

#Accuracy factor
af <- exp(sqrt(af_sum/nrow(shelf_life_data_af)))

## ---------------------------Discrepancy---------------------------------------

percent_discrep <- (af - 1)*100
percent_discrep

## ---------------------------Bias----------------------------------------------

bf_sum <- 0

for(i in 1:nrow(shelf_life_data_af)){
  d <- shelf_life_data_af$sample_id[i]
  add <- (med_shelf_life_af$median_ln[med_shelf_life_af$day == d] -
            shelf_life_data_af$geom_ln_conc[shelf_life_data_af$sample_id == d])
  bf_sum <- bf_sum + add
}

#Bias factor
bf <- exp(bf_sum/nrow(shelf_life_data_af))

## ---------------------------Percent Bias--------------------------------------

percent_bias <- (-1)*(exp(abs(log(bf))) - 1)*100
percent_bias

## ---------------------------Parameter values for packages---------------------

d18_summary <- packages_summary %>%
  filter(day == "18") %>%
  mutate(abov_thresh = ifelse(total_count > 8.0, "y", "n"))

d18_param_summary <- d18_summary %>%
  group_by(abov_thresh) %>%
  summarize(n = n()
,            mean_init_count = round(mean(initial_count_apc), 2),
            mean_temp = round(mean(storage_temp), 2),
            mean_strains = round(mean(tot_strains), 2),
            mean_nmax = round(mean(final_nmax), 2),
            mean_mmax = round(mean(mean_mumax), 2),
            mean_lg = round(mean(mean_lag), 2),
            prop_s120116 = round(sum(s12.0116)/n(), 2),
            prop_s120132 = round(sum(s12.0132)/n(), 2),
            prop_s120141 = round(sum(s12.0141)/n(), 2),
            prop_s120166 = round(sum(s12.0166)/n(), 2),
            prop_s120180 = round(sum(s12.0180)/n(), 2),
            prop_s120184 = round(sum(s12.0184)/n(), 2))


