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

strain_param <- read.csv("model_inputs/primary_model_parameters_strain_averaged_2023_02_07.csv", header = TRUE)

mumax_secon <- read.csv("model_inputs/secondary_model_mumax_strain_averaged_2023_02_08.csv", header = TRUE)
lag_secon <- read.csv("model_inputs/secondary_model_lag_strain_averaged_2023_02_08.csv", header = TRUE)

initial_count <- read.csv("model_inputs/initial_distribution_2023_03_10.csv", header = TRUE)
st_frequency <- read.csv("model_inputs/st_frequency_2022_05_18.csv", header = TRUE)

#Change strain parameters to lowercase 
strain_param$isolate <- tolower(strain_param$isolate)

mumax_secon$isolate <- tolower(mumax_secon$isolate)
lag_secon$isolate <- tolower(lag_secon$isolate)

#Initial Count Distribution: Log transforming initial count, and fitting a uniform distribution to the data
initial_count_dist <- function(n){
  ans <- runif(n, min = initial_count$min, max = initial_count$max)
  return(ans)
}

#ST Distribution
isolates <- st_frequency$cida_isolates
st_prev <- st_frequency$frequency

st_dist <- function(n) {
ans <- sample(isolates, size = n, prob = st_prev, replace = F)
return(ans)
}

#Storage Temperature Distribution
temp_location <- 4
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

#pop_dist <- function(n){
#  ans <- runif(n, min = 1, max = 6)
#  ans_rounded <- round(ans, digits = 0)
#  return(ans_rounded)
#}

#Primary Growth Models 
#Gompertz
gompertz <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10n0 + (log10nmax - log10n0) * exp(-exp(mumax * exp(1) * (lag - day) / ((log10nmax - log10n0) * log(10)) + 1))
}

#Secondary model 

#Mumax
NewMu <- function(temp, b, temp_min) {
  sqrt_mu <- b*(temp - temp_min)
  ans <- sqrt_mu^2
  return(ans)
}

#Lag phase
NewLag <- function(temp, b, temp_min) {
  sqrt_lag <- b*(temp - temp_min)
  ans <- sqrt_lag^2
  return(ans)
}

## ---------------------------Simulation Overview-------------------------------

n_lots <- 1
n_packages <- 1000
n_sim <- n_packages*n_lots
n_day <- 22

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, number of strains, STs and temperature to each package
packages_sim <- data.frame(matrix(, nrow = n_sim, ncol = 4))
colnames(packages_sim) <- c("package", "initial_count_apc", "storage_temp", "tot_strains")

packages_sim$package <- 1:n_packages
packages_sim$initial_count_apc <- initial_count_dist(n_sim)
packages_sim$storage_temp <- temp_dist(n_sim)
packages_sim$tot_strains <- 3 #pop_dist(n_sim)

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
  strain_index <- which(packages_sim_exp$isolate[i] == lag_secon$isolate)
  strain_lag <- NewLag(packages_sim_exp$storage_temp[i], lag_secon$b_lag[strain_index], lag_secon$temp_min_lag[strain_index])
  packages_sim_exp$lag[i] <- strain_lag
}

#Assign mumax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == mumax_secon$isolate)
  strain_mu <- NewMu(packages_sim_exp$storage_temp[i], mumax_secon$b_mu[strain_index], mumax_secon$temp_min_mu[strain_index])
  packages_sim_exp$mumax[i] <- strain_mu
}

#Assign nmax to each package
for(i in 1:nrow(packages_sim_exp)){
  strain_index <- which(packages_sim_exp$isolate[i] == strain_param$isolate & strain_param$temp == 6)
  strain_nmax <- strain_param$nmax[strain_index]
  packages_sim_exp$nmax[i] <- strain_nmax
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
  summarize(spoiled = sum(total_count > 7.6)/n_packages)

h18 <- hist(day18_sim_output$total_count, 
            breaks = 15,
            xlab = "Bacterial concentration, log10(CFU/g)", 
            ylab = "Density", main = "New model output on day 18 of shelf life", 
            prob = TRUE, col = "white")
abline(v = 7.6, lty = 4)
legend(
  'topleft', ncol = 1,
  legend = c(
    'Percent of packages spoiled',  (d18_spoiled$spoiled*100)))

mean_count_summary <- shelf_life_sim %>%
  group_by(day, package) %>%
  summarize(total_count = log10(sum(10^(count)))) %>%
  group_by(day) %>%
  summarize(mean_total_count = mean(total_count), sd_total_count = sd(total_count))


## --------------------------------Threshold------------------------------------

sa <- shelf_life_sim %>%
  group_by(package,day) %>%
  summarize(total_count = log10(sum(10^(count)))) %>%
  group_by(day) %>%
  summarize(prop_spoiled = sum(total_count > 7.6)/n_packages,
            variance = var(total_count)) %>%
  mutate(adj = 3) %>%
  mutate(param = "no_of_strains")

## --------------------------------Export---------------------------------------

date <- Sys.Date()
date <- gsub("-", "_", date)
file_name <- paste0("sensitivity_analysis/strain_no/output/sa_strain_no_", date,".csv")

write.table(sa, file_name, sep = ",", col.names = !file.exists(file_name), append = T, row.names = FALSE)
