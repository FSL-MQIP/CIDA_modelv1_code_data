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

strain_param <- read.csv("primary_model_parameters_strain_by_rep_2022_08_26.csv", header = TRUE)
apc_param <- read.csv("primary_model_parameters_apc_by_rep_2022_09_01.csv", header = TRUE)

mumax_secon <- read.csv("secondary_model_mumax_strain_by_rep_2022-09-01.csv", header = TRUE)
lag_secon <- read.csv("secondary_model_lag_strain_by_rep_2022_09_01.csv", header = TRUE)
mumax_secon$b_mu <- 1.2*mumax_secon$b_mu
lag_apc_secon <- read.csv("secondary_model_lag_apc_by_rep2022_09_01.csv", header = TRUE)

initial_count <- read.csv("initial_distribution_2022_08_24.csv", header = TRUE)
st_frequency <- read.csv("st_frequency_2022_05_18.csv", header = TRUE)

#Change strain parameters to lowercase 
strain_param$isolate <- tolower(strain_param$isolate)
apc_param$isolate <- tolower(apc_param$isolate)

mumax_secon$isolate <- tolower(mumax_secon$isolate)
lag_secon$isolate <- tolower(lag_secon$isolate)

#Initial Count Distribution: Log transforming initial count, and fitting a normal distribution to the data
initial_count_dist <- function(n){
  ans <- runif(n, min = initial_count$min, max = initial_count$max)
  return(ans)
}

#ST Distribution
isolates <- st_frequency$cida_isolates
st_prev <- st_frequency$frequency

st_dist <- function(n) {
ans <- sample(isolates, size = n, prob = st_prev, replace = T)
return(ans)
}

#Storage Temperature Distribution
temp_location <- 4
temp_shape <- 2.31

temp_dist <- function(n) {
  ans <- rlaplace(n, temp_location, temp_shape)
  for(i in 1:n)
  while (ans[i] < 0 | ans[i] > 10) {
    ans[i] <- rlaplace(1, temp_location, temp_shape)
  }
  return(ans)
}

#Primary Growth Models 


#Baranyi
baranyi <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10nmax + log10((-1 + exp(mumax * lag) + exp(mumax * day))/(exp(mumax * day) - 1 + exp(mumax * lag) * 10^(log10nmax - log10n0)))
}

#Baranyi without lag
baranyi_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10nmax - log10(1 + (10^(log10nmax - log10n0) - 1) * exp(-mumax * day))
}

#Buchanan
buchanan <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10n0 + (day >= lag) * (day <= (lag + (log10nmax - log10n0) * log(10) / mumax)) * mumax * (day - lag) / log(10) + (day >= lag) * (day>(lag + (log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

#Buchanan without lag
buchanan_no_lag <- function(day, log10n0, log10nmax, mumax){
  log10n <- log10n0 + (day <= ((log10nmax - log10n0) * log(10) / mumax)) * mumax * day / log(10) + (day > ((log10nmax - log10n0) * log(10) / mumax)) * (log10nmax - log10n0)
}

#Gompertz
gompertz <- function(day, log10n0, log10nmax, mumax, lag){
  log10n <- log10n0 + (log10nmax - log10n0) * exp(-exp(mumax * exp(1) * (lag - day) / ((log10nmax - log10n0) * log(10)) + 1))
}


log10N_with_lag <- function(day, log10n0, log10nmax, mumax, lag, model) {
  ans <- numeric(length = length(day))
  for(i in seq_along(day)){
    if (model[i] == "baranyi") {
      ans[i] <- baranyi(day[i], log10n0[i], log10nmax[i], mumax[i], lag[i])
    }
    else if(model[i] == 'buchanan') {
      ans[i] <- buchanan(day[i], log10n0[i], log10nmax[i], mumax[i], lag[i])
    }
    else if(model[i] == "gompertz"){
      ans[i] <- gompertz(day[i], log10n0[i], log10nmax[i], mumax[i], lag[i])
    } 
    else{
      stop(paste0(model, " is not a valid model name. Must be one of Buchanan, Baranyi, Gompertz."))
    }
  }
  return(ans)
  
}

log10N_without_lag <- function(day, log10n0, log10nmax, mumax, model) {
  ans <- numeric(length = length(day))
  for(i in seq_along(day)){
    if (model[i] == "baranyi_no_lag") {
      ans[i] <- baranyi_no_lag(day[i], log10n0[i], log10nmax[i], mumax[i])
    }
    else if(model[i] == 'buchanan_no_lag') {
      ans[i] <- buchanan_no_lag(day[i], log10n0[i], log10nmax[i], mumax[i])
    }
    else{
      stop(paste0(model, " is not a valid model name. Must be one of Buchanan, Baranyi, Gompertz."))
    }
  }
  return(ans)
  
}

#Secondary model 

NewMu <- function(temp, b, temp_min) {
  sqrt_mu <- b*(temp - temp_min)
  ans <- sqrt_mu^2
  return(ans)
}

NewLag <- function(temp, b, temp_min) {
  sqrt_lag <- b*(temp - temp_min)
  ans <- sqrt_lag^2
  return(ans)
}

## ---------------------------Simulation Overview-------------------------------

n_packages <- 10000
n_sim <- n_packages
n_day <- 22

## ---------------Assigning Model Inputs----------------------------------------

#Assigning initial count, ST and temperature to each package
packages_sim <- data.frame(matrix(, nrow = n_sim, ncol = 9))
colnames(packages_sim) <- c("package", "initial_count", "storage_temp", "isolate", "isolate_rep", "model", "lag", "mu", "nmax")

packages_sim$package <- rep(1:n_packages)
packages_sim$initial_count <- initial_count_dist(n_sim)
packages_sim$storage_temp <- temp_dist(n_sim)
packages_sim$isolate <- st_dist(n_sim)
packages_sim$isolate_rep <- sample(1:3, n_sim, replace = TRUE)

#Assigning a growth model to each package. The selected model will be the best fitting model at 6˚C 
for(i in 1:n_sim){
  strain_index <- which(packages_sim$isolate[i] == strain_param$isolate & packages_sim$isolate_rep[i] == strain_param$rep & 6 == strain_param$temp)
  packages_sim$model[i] <- strain_param$model_type[strain_index]
}

#Assign lag to each package

for(i in 1:n_sim){
  rep <- packages_sim$isolate_rep[i]
  if(packages_sim$isolate[i] == "s12-0180"){
    rep <- sample(c(1, 3), 1, replace = TRUE)
  } else if (packages_sim$isolate[i] == "s12-0184"){
    rep <- 1
  } else if (packages_sim$isolate[i] == "s12-0132"){
    rep <- 3
  }  else {
    rep <- rep
  }
  apc_lag <- 0
  strain_index <- which(packages_sim$isolate[i] == lag_secon$isolate & rep == lag_secon$rep)
  strain_lag <- NewLag(packages_sim$storage_temp[i], lag_secon$b_mu[strain_index], lag_secon$temp_min_mu[strain_index])
  strain_lag <- ifelse(packages_sim$isolate[i] == "s12-0166", 0, strain_lag)
  packages_sim$lag[i] <- runif(n = 1, min = min(strain_lag, apc_lag), max = max(strain_lag, apc_lag))
  packages_sim$lag[i] <- packages_sim$lag[i]
}

packages_sim$lag <- ifelse(packages_sim$model %in% c("baranyi_no_lag", "buchanan_no_lag"), NA, packages_sim$lag)

#Assign mumax to each package
for(i in 1:n_sim){
  rep <- ifelse(packages_sim$isolate[i] == "s12-0141", 1, packages_sim$isolate_rep[i])
  strain_index <- which(packages_sim$isolate[i] == mumax_secon$isolate & rep == mumax_secon$rep)
  package_mu <- NewMu(packages_sim$storage_temp[i], mumax_secon$b_mu[strain_index], mumax_secon$temp_min_mu[strain_index])
  packages_sim$mu[i] <- package_mu
}

#Assign nmax to each package
for(i in 1:n_sim){
  rep <- packages_sim$isolate_rep[i]
strain_index <- which(packages_sim$isolate[i] == strain_param$isolate & rep == strain_param$rep & 6 == strain_param$temp)
apc_index <- which(packages_sim$isolate[i] == apc_param$isolate & rep == apc_param$rep & 6 == apc_param$temp)
strain_nmax <- strain_param$nmax[strain_index]
apc_nmax <- apc_param$nmax[strain_index]
packages_sim$nmax[i] <- runif(n = 1, strain_nmax, apc_nmax)
}

## ----------------Assessing Growth Over Shelf Life-----------------------------

models_with_lag <- c("baranyi", "buchanan", "gompertz")

shelf_life_sim <- packages_sim[rep(seq_len(nrow(packages_sim)), each = 22), ]
rownames(shelf_life_sim) <- 1:nrow(shelf_life_sim)
shelf_life_sim["day"] <- rep(1:n_day, times = n_sim)
shelf_life_sim["count"] <- NA

shelf_life_sim[is.na(shelf_life_sim$lag), ] <- shelf_life_sim[is.na(shelf_life_sim$lag), ] %>% 
  mutate(count = log10N_without_lag(day, initial_count, nmax, mu, model))

shelf_life_sim[!is.na(shelf_life_sim$lag), ] <- shelf_life_sim[!is.na(shelf_life_sim$lag), ] %>% 
  mutate(count = log10N_with_lag(day, initial_count, nmax, mu, lag, model))

## --------------------------------Plotting Outcome-----------------------------

day17_sim_output <- shelf_life_sim %>%
  filter(day == 17)

h17 <- hist(day17_sim_output$count)
h17


day22_sim_output <- shelf_life_sim %>%
  filter(day == 22)

h22 <- hist(day22_sim_output$count)
h22

## --------------------------------Threshold------------------------------------

sa_temp <- shelf_life_sim %>%
  group_by(day) %>%
  summarise(prop_fail = sum(count > 7.6)/n_packages) %>%
  mutate(adj = 1.2) %>%
  mutate(param = "b_mumax")

## --------------------------------Export---------------------------------------

date <- Sys.Date()
date <- gsub("-", "_", date)
file_name <- paste0("output/sa_b_mumax_", date,".csv")
  
write.table(sa_temp, file_name, sep = ",", col.names = !file.exists(file_name), append = T, row.names = FALSE)


