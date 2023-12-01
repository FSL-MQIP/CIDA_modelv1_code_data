## ---------------------------Title---------------------------------------------
# CIDA Shelf Life Model

## ---------------------------Description---------------------------------------
# Project: CIDA Spinach 

# Description: 

## ---------------------------Packages------------------------------------------

#   Loading packages
library(dplyr); library(MASS); library(VGAM)

## ---------------------------Seed----------------------------------------------

#Set seed
set.seed(1)

## ---------------------------Inputs--------------------------------------------

#Load in data

strain_param <- read.csv("model_inputs/primary_model_parameters_strain_averaged.csv", header = TRUE)

mumax_secon <- read.csv("model_inputs/secondary_model_mumax_strain_averaged.csv", header = TRUE)

st_frequency <- read.csv("model_inputs/st_frequency.csv", header = TRUE)

initial_count <- read.csv("model_inputs/initial_distribution.csv", header = TRUE)

#Change strain parameters to lowercase 
strain_param$isolate <- tolower(strain_param$isolate)

mumax_secon$isolate <- tolower(mumax_secon$isolate)

#Change the formatting of a variable in the ST frequency 
st_frequency$growth_isolates <- substr(st_frequency$growth_isolates, 1, 8)
st_frequency$growth_isolates <- tolower(st_frequency$growth_isolates)

#Initial Count Distribution: Log transforming initial count, and fitting a uniform distribution to the data
initial_count_dist <- function(n){
  ans <- runif(n, min = initial_count$min, max = initial_count$max)
  return(ans)
}

#ST Distribution
isolates <- st_frequency$growth_isolates
st_prev <- st_frequency$frequency

st_dist <- function(n) {
ans <- sample(isolates, size = n, prob = st_prev, replace = F)
return(ans)
}

#Storage Temperature Distribution
temp_location <- 4.06
temp_shape <- 2.31

temp_dist <- function(n) {
  ans <- rlaplace(n, temp_location, temp_shape)
  for(i in 1:n)
  while (ans[i] < 2 | ans[i] > 14) {
    ans[i] <- rlaplace(1, temp_location, temp_shape)
  }
  return(ans)
}

#Number of strains per package

pop_dist <- function(n){
  ans <- runif(n, min = 1, max = 6)
  ans_rounded <- round(ans, digits = 0)
  return(ans_rounded)
}

#Primary Growth Models 
#Gompertz
gompertz <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10n0 + (log10n0 < log10nmax)*(log10nmax - log10n0) * exp(-exp(mumax * exp(1) * (lag - day) / ((log10nmax - log10n0) * log(10)) + 1))
}

#Secondary model 

#Mumax
NewMu_general <- function(temp, b, temp_min){
    if(temp > temp_min){
      sqrt_mu <- b*(temp - temp_min)
      ans <- sqrt_mu^2
      return(ans)
    } else {
    ans <- 0 
    return(ans)
  }
}

NewMu_s12_0141 <- function(temp, b, temp_min){
  if(temp < temp_min){
    sqrt_mu <- b*(temp - temp_min)
    ans <- sqrt_mu^2
    return(ans)
  } else {
    ans <- 0 
    return(ans)
  }
}

NewMu <- function(temp, b, temp_min, isolate){
  if(isolate != "s12-0141"){
    return(NewMu_general(temp, b, temp_min))
  } else {
    return(NewMu_s12_0141(temp, b, temp_min))
  }
}

lag_tmin_pre <- mumax_secon$temp_min_mu[mumax_secon$isolate == "s12-0116"]
lag_tmin <- lag_tmin_pre + 60*(lag_tmin_pre/100)

NewLag <- function(temp, lag) {
  ans <- lag*((6 - lag_tmin)/(temp - lag_tmin))^2
  return(ans)
}

## ---------------------------Simulation Overview-------------------------------

n_lots <- 1
n_packages <- 1000
n_sim <- n_packages*n_lots
n_day <- 22

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, number of strains, STs and temperature to each package
packages_sim <- data.frame(matrix(nrow = n_sim, ncol = 4))
colnames(packages_sim) <- c("package", "initial_count_apc", "storage_temp", "tot_strains")

packages_sim$package <- 1:n_packages
packages_sim$initial_count_apc <- initial_count_dist(n_sim)
packages_sim$storage_temp <- temp_dist(n_sim)
packages_sim$tot_strains <- pop_dist(n_sim)

#Expanding the dataframe to store observations for each strain and day of observation
packages_sim_exp <- packages_sim %>%
  group_by(package) %>%
  slice(rep(1, times = tot_strains)) %>%
  mutate(isolate = NA, initial_count = NA, lag = NA, mumax = NA, nmax = NA) 

#Assigning isolates to each package, with replacement
for(i in 1:n_packages){
  index <- which(packages_sim_exp$package == i)
  strain_per_package <- unique(packages_sim_exp$tot_strains[index])
  packages_sim_exp$isolate[index] <- st_dist(strain_per_package) 
  packages_sim_exp$initial_count[index] <- rep(log10((10^(packages_sim$initial_count_apc[i]))/strain_per_package), times = strain_per_package)
}

#Assign lag to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == strain_param$isolate & 6 == strain_param$temp)
  strain_lag <- NewLag(packages_sim_exp$storage_temp[i], strain_param$lag[strain_index])
  packages_sim_exp$lag[i] <- strain_lag
}

#Assign mumax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == mumax_secon$isolate)
  strain_mu <- NewMu(packages_sim_exp$storage_temp[i], mumax_secon$b_mu[strain_index], mumax_secon$temp_min_mu[strain_index], packages_sim_exp$isolate[i])
  packages_sim_exp$mumax[i] <- strain_mu
}

#Assign nmax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == strain_param$isolate & strain_param$temp == 6)
  strain_nmax <- strain_param$nmax[strain_index]
  packages_sim_exp$nmax[i] <- strain_nmax
}

#Calculate the total nmax of a package 
packages_sim_exp$final_nmax <- vector(mode = "integer", length = nrow(packages_sim_exp))

for(i in 1:n_packages){
  index <- which(packages_sim_exp$package == i)
  packages_nmax <- packages_sim_exp$nmax[index]
  final_nmax <- log10(sum(10^packages_nmax))
  packages_sim_exp$final_nmax[index] <- final_nmax
}

## ----------------Assessing Growth Over Shelf Life-----------------------------

#Creating a dataframe, to assess microbial growth over shelf life 
shelf_life_sim <- packages_sim_exp %>%
  group_by(package, isolate) %>%
  slice(rep(1, times = 22)) %>%
  mutate(day = 1:22)
  
shelf_life_sim <- shelf_life_sim %>%
  mutate(count = gompertz(day, initial_count, nmax, mumax, lag))

## --------------------Plotting & Analyzing Outcome-----------------------------

day18_sim_output <- shelf_life_sim %>%
  filter(day == 18) %>%
  group_by(package) %>%
  summarize(total_count = log10(sum(10^(count))))

d18_spoiled <- day18_sim_output %>%
  summarize(spoiled = sum(total_count > 8.0)/n_packages)

h18 <- hist(day18_sim_output$total_count, 
            breaks = 15,
            xlab = "Bacterial concentration, log10(CFU/g)", 
            ylab = "Density", main = "Simulation without secondary model for lag", 
            prob = TRUE, col = "white")
abline(v = 8.0, lty = 4)
legend(
  'topleft', ncol = 1,
  legend = c(
    'Percent of packages spoiled on Day 18',  (d18_spoiled$spoiled*100)))

percent_packages_nmax <- shelf_life_sim %>%
  group_by(day, package) %>%
  mutate(total_count = log10(sum(10^(count)))) %>%
  slice(1) %>%
  mutate(at_nmax = round(final_nmax, 1)  == round(total_count, 1)) %>%
  group_by(day) %>%
  summarize(percent_packs_nmax = (sum(at_nmax == TRUE)/n_packages)*100)

mode <- function(x) {
  t <- table(x)
  names(t)[ which.max(t) ]
}

count_summary_by_day <- shelf_life_sim %>%
  group_by(day, package) %>%
  summarize(total_count = log10(sum(10^(count)))) %>%
  mutate(rounded_total_count = round(total_count, 1)) %>%
  group_by(day) %>%
  summarize(mean = mean(total_count), 
            sd_total_count = sd(total_count), 
            mode = mode(rounded_total_count),
            median = median(total_count))

counts_summary_by_package <- shelf_life_sim %>%
  group_by(day, package) %>%
  summarize(total_count = log10(sum(10^(count))))

## --------------------------------Threshold------------------------------------

sa_1 <- shelf_life_sim %>%
  group_by(package,day) %>%
  summarize(total_count = log10(sum(10^(count)))) %>%
  group_by(day) %>%
  summarize(prop_spoiled = sum(total_count > 8.0)/n_packages) %>%
  mutate(adj = 0.4) %>%
  mutate(param = "lag_tmin")

sa_2 <- shelf_life_sim %>%
  group_by(day, package) %>%
  summarize(total_count = log10(sum(10^(count)))) %>%
  mutate(rounded_total_count = round(total_count, 1)) %>%
  group_by(day) %>%
  summarize(mean = mean(total_count), 
            sd_total_count = sd(total_count), 
            mode = mode(rounded_total_count),
            median = median(total_count)) %>%
  mutate(adj = 0.4) %>%
  mutate(param = "lag_tmin")

## --------------------------------Export---------------------------------------

file_name_1 <- "sensitivity_analysis/tmin_lag/output/sa_lag_tmin_threshold.csv"
file_name_2 <- "sensitivity_analysis/tmin_lag/output/sa_lag_tmin_count.csv"

write.table(sa_1, file_name_1, sep = ",", col.names = !file.exists(file_name_1), append = T, row.names = FALSE)
write.table(sa_2, file_name_2, sep = ",", col.names = !file.exists(file_name_2), append = T, row.names = FALSE)

