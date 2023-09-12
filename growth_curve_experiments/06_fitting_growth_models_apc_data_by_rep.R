##  ----------------------------Title-------------------------------------------
#   Fitting models to the SM-ACPC growth data

##  --------------------------Description---------------------------------------
#   Project: CIDA Spinach 

#   Script description: Primary growth models, Baranyi (with and without lag), Buchanan (with and without lag) and Gompertz) will be fit to the SM-APC growth data. 

##  --------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(minpack.lm)

##  --------------------------Data----------------------------------------------
#Read in data 
strain_apc_data <- read.csv("data/wrangled/apc_strain_averaged_wrangled_03.csv", header = TRUE)
param_estimates <- read.csv("outputs/parameters/initial_parameter_estimates_apc_by_rep.csv", header = TRUE)

## ------------------Defining Primary Models------------------------------------
# Defining the primary growth functions. The Baranyi (with and without lag), Buchanan (with and without lag), and Gompertz model are found in the nlsMicrobio package. 
#Link to nlsMicrobio on Cran: https://cran.r-project.org/web/packages/nlsMicrobio/index.html. Mumax is in log (base e).
#The gompertz model here is the modified Gompertz model (Zwietering et al., 1990)

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

#End of defining primary models

## ---------------------------Data Wrangling------------------------------------
#Data is split into apc growth data and APC growth data 
apc_data <- strain_apc_data %>%
  dplyr::filter(media == "BHI")

#Subset apc data by isolate and rep
b116_6_1 <- filter(apc_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B2")
b116_6_2 <- filter(apc_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B3")
b116_6_3 <- filter(apc_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B6")

b132_6_1 <- filter(apc_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B2")
b132_6_2 <- filter(apc_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B3")
b132_6_3 <- filter(apc_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B6")

b141_6_1 <- filter(apc_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B2")
b141_6_2 <- filter(apc_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B3")
b141_6_3 <- filter(apc_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B6")

b166_6_1 <- filter(apc_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B4")
b166_6_2 <- filter(apc_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B5")
b166_6_3 <- filter(apc_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B7")

b180_6_1 <- filter(apc_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B4")
b180_6_2 <- filter(apc_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B5")
b180_6_3 <- filter(apc_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B7")

b184_6_1 <- filter(apc_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B4")
b184_6_2 <- filter(apc_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B5")
b184_6_3 <- filter(apc_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B7")

b116_10 <- filter(apc_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "10")
b132_10 <- filter(apc_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "10")
b141_10 <- filter(apc_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "10")
b166_10 <- filter(apc_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "10")
b180_10 <- filter(apc_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "10")
b184_10 <- filter(apc_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "10")

#Making a dataframe to store the fitted model parameter estimates for the isolates 
no_of_models <- 5
fitted_param <- matrix(nrow = (24 * no_of_models), ncol = 13)
fitted_param <- data.frame(fitted_param)
colnames(fitted_param) <- c("isolate", "pop", "model_type", "n0", "lag", "mumax", "nmax", "temp", "aic", "bic", "fit", "method", "rep")

#End of data wrangling

##---------------------------S12-0116 Model Fitting, 6C-------------------------
#Rep 1
b116_6_1_index <- which(param_estimates$isolate == "S12-0116" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b116_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b116_6_1, 
                               start=list(
                                 log10n0 = param_estimates$n0[b116_6_1_index], 
                                 log10nmax = param_estimates$nmax[b116_6_1_index], 
                                 mumax = (param_estimates$mumax[b116_6_1_index]*2.303), 
                                 lag = param_estimates$lag[b116_6_1_index]), 
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[1] <-"S12-0116"
fitted_param$pop[1] <-"apc"
fitted_param$model_type[1] <-"baranyi"
fitted_param$n0[1] <- summary(b116_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[1] <- summary(b116_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[1] <- summary(b116_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[1] <- summary(b116_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[1] <- "6"
fitted_param$aic[1] <- AIC(b116_baranyi_nls_lm_6_1)
fitted_param$bic[1] <- BIC(b116_baranyi_nls_lm_6_1)
fitted_param$fit[1] <- "yes"
fitted_param$method[1] <- "nls_lm"
fitted_param$rep[1] <- "1"

#baranyi no lag
b116_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b116_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b116_6_1_index], 
                                      log10nmax = param_estimates$nmax[b116_6_1_index], 
                                      mumax = (param_estimates$mumax[b116_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[2] <-"S12-0116"
fitted_param$pop[2] <-"apc"
fitted_param$model_type[2] <-"baranyi_no_lag"
fitted_param$n0[2] <- summary(b116_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[2] <- summary(b116_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[2] <- summary(b116_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[2] <- "6"
fitted_param$aic[2] <- AIC(b116_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[2] <- BIC(b116_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[2] <- "yes"
fitted_param$method[2] <- "nls_lm"
fitted_param$rep[2] <- "1"

#buchanan
b116_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b116_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b116_6_1_index], 
                                log10nmax = param_estimates$nmax[b116_6_1_index], 
                                mumax = (param_estimates$mumax[b116_6_1_index]*2.303), 
                                lag = param_estimates$lag[b116_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[3] <-"S12-0116"
fitted_param$pop[3] <-"apc"
fitted_param$model_type[3] <-"buchanan"
fitted_param$n0[3] <- summary(b116_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[3] <- summary(b116_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[3] <- summary(b116_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[3] <- summary(b116_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[3] <- "6"
fitted_param$aic[3] <- AIC(b116_buchanan_nls_lm_6_1)
fitted_param$bic[3] <- BIC(b116_buchanan_nls_lm_6_1)
fitted_param$fit[3] <- "yes"
fitted_param$method[3] <- "nls_lm"
fitted_param$rep[3] <- "1"

#buchanan no lag
b116_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b116_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b116_6_1_index], 
                                       log10nmax = param_estimates$nmax[b116_6_1_index], 
                                       mumax = (param_estimates$mumax[b116_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[4] <-"S12-0116"
fitted_param$pop[4] <-"apc"
fitted_param$model_type[4] <-"buchanan_no_lag"
fitted_param$n0[4] <- summary(b116_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[4] <- summary(b116_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[4] <- summary(b116_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[4] <- "6"
fitted_param$aic[4] <- AIC(b116_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[4] <- BIC(b116_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[4] <- "yes"
fitted_param$method[4] <- "nls_lm"
fitted_param$rep[4] <- "1"

#gompertz
b116_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b116_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b116_6_1_index], 
                                log10nmax = param_estimates$nmax[b116_6_1_index], 
                                mumax = (param_estimates$mumax[b116_6_1_index]*2.303),
                                lag = param_estimates$lag[b116_6_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[5] <-"S12-0116"
fitted_param$pop[5] <-"apc"
fitted_param$model_type[5] <-"gompertz"
fitted_param$n0[5] <- summary(b116_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[5] <- summary(b116_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[5] <- summary(b116_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[5] <- summary(b116_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[5] <- "6"
fitted_param$aic[5] <- AIC(b116_gompertz_nls_lm_6_1)
fitted_param$bic[5] <- BIC(b116_gompertz_nls_lm_6_1)
fitted_param$fit[5] <- "yes"
fitted_param$method[5] <- "nls_lm"
fitted_param$rep[5] <- "1"

#Rep 2
b116_6_2_index <- which(param_estimates$isolate == "S12-0116" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b116_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b116_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b116_6_2_index], 
                                 log10nmax = param_estimates$nmax[b116_6_2_index], 
                                 mumax = (param_estimates$mumax[b116_6_2_index]*2.303), 
                                 lag = param_estimates$lag[b116_6_2_index]), 
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[6] <-"S12-0116"
fitted_param$pop[6] <-"apc"
fitted_param$model_type[6] <-"baranyi"
fitted_param$n0[6] <- summary(b116_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[6] <- summary(b116_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[6] <- summary(b116_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[6] <- summary(b116_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[6] <- "6"
fitted_param$aic[6] <- AIC(b116_baranyi_nls_lm_6_2)
fitted_param$bic[6] <- BIC(b116_baranyi_nls_lm_6_2)
fitted_param$fit[6] <- "yes"
fitted_param$method[6] <- "nls_lm"
fitted_param$rep[6] <- "2"

#baranyi no lag
b116_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b116_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b116_6_2_index], 
                                        log10nmax = param_estimates$nmax[b116_6_2_index], 
                                        mumax = (param_estimates$mumax[b116_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[7] <-"S12-0116"
fitted_param$pop[7] <-"apc"
fitted_param$model_type[7] <-"baranyi_no_lag"
fitted_param$n0[7] <- summary(b116_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[7] <- summary(b116_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[7] <- summary(b116_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[7] <- "6"
fitted_param$aic[7] <- AIC(b116_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[7] <- BIC(b116_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[7] <- "yes"
fitted_param$method[7] <- "nls_lm"
fitted_param$rep[7] <- "2"

#buchanan
b116_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b116_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b116_6_2_index], 
                                  log10nmax = param_estimates$nmax[b116_6_2_index], 
                                  mumax = (param_estimates$mumax[b116_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b116_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[8] <-"S12-0116"
fitted_param$pop[8] <-"apc"
fitted_param$model_type[8] <-"buchanan"
fitted_param$n0[8] <- summary(b116_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[8] <- summary(b116_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[8] <- summary(b116_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[8] <- summary(b116_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[8] <- "6"
fitted_param$aic[8] <- AIC(b116_buchanan_nls_lm_6_2)
fitted_param$bic[8] <- BIC(b116_buchanan_nls_lm_6_2)
fitted_param$fit[8] <- "yes"
fitted_param$method[8] <- "nls_lm"
fitted_param$rep[8] <- "2"

#buchanan no lag
b116_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b116_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b116_6_2_index], 
                                         log10nmax = param_estimates$nmax[b116_6_2_index], 
                                         mumax = (param_estimates$mumax[b116_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[9] <-"S12-0116"
fitted_param$pop[9] <-"apc"
fitted_param$model_type[9] <-"buchanan_no_lag"
fitted_param$n0[9] <- summary(b116_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[9] <- summary(b116_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[9] <- summary(b116_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[9] <- "6"
fitted_param$aic[9] <- AIC(b116_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[9] <- BIC(b116_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[9] <- "yes"
fitted_param$method[9] <- "nls_lm"
fitted_param$rep[9] <- "2"

#gompertz
b116_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b116_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b116_6_2_index], 
                                  log10nmax = param_estimates$nmax[b116_6_2_index], 
                                  mumax = (param_estimates$mumax[b116_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b116_6_2_index]), 
                                lower = c(0, 0, 0, 0), 
                                control = nls.lm.control(maxiter = 150))

fitted_param$isolate[10] <-"S12-0116"
fitted_param$pop[10] <-"apc"
fitted_param$model_type[10] <-"gompertz"
fitted_param$n0[10] <- summary(b116_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[10] <- summary(b116_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[10] <- summary(b116_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[10] <- summary(b116_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[10] <- "6"
fitted_param$aic[10] <- AIC(b116_gompertz_nls_lm_6_2)
fitted_param$bic[10] <- BIC(b116_gompertz_nls_lm_6_2)
fitted_param$fit[10] <- "yes_maxiter_error"
fitted_param$method[10] <- "nls_lm"
fitted_param$rep[10] <- "2"

#Rep 3
b116_6_3_index <- which(param_estimates$isolate == "S12-0116" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b116_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b116_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b116_6_3_index], 
                                 log10nmax = param_estimates$nmax[b116_6_3_index], 
                                 mumax = (param_estimates$mumax[b116_6_3_index]*2.303), 
                                 lag = param_estimates$lag[b116_6_3_index]), 
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[11] <-"S12-0116"
fitted_param$pop[11] <-"apc"
fitted_param$model_type[11] <-"baranyi"
fitted_param$n0[11] <- summary(b116_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[11] <- summary(b116_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[11] <- summary(b116_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[11] <- summary(b116_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[11] <- "6"
fitted_param$aic[11] <- AIC(b116_baranyi_nls_lm_6_3)
fitted_param$bic[11] <- BIC(b116_baranyi_nls_lm_6_3)
fitted_param$fit[11] <- "yes"
fitted_param$method[11] <- "nls_lm"
fitted_param$rep[11] <- "3"

#baranyi no lag
b116_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b116_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b116_6_3_index], 
                                        log10nmax = param_estimates$nmax[b116_6_3_index], 
                                        mumax = (param_estimates$mumax[b116_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[12] <-"S12-0116"
fitted_param$pop[12] <-"apc"
fitted_param$model_type[12] <-"baranyi_no_lag"
fitted_param$n0[12] <- summary(b116_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[12] <- summary(b116_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[12] <- summary(b116_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[12] <- "6"
fitted_param$aic[12] <- AIC(b116_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[12] <- BIC(b116_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[12] <- "yes"
fitted_param$method[12] <- "nls_lm"
fitted_param$rep[12] <- "3"

#buchanan
b116_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b116_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b116_6_3_index], 
                                  log10nmax = param_estimates$nmax[b116_6_3_index], 
                                  mumax = (param_estimates$mumax[b116_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b116_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[13] <-"S12-0116"
fitted_param$pop[13] <-"apc"
fitted_param$model_type[13] <-"buchanan"
fitted_param$n0[13] <- summary(b116_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[13] <- summary(b116_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[13] <- summary(b116_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[13] <- summary(b116_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[13] <- "6"
fitted_param$aic[13] <- AIC(b116_buchanan_nls_lm_6_3)
fitted_param$bic[13] <- BIC(b116_buchanan_nls_lm_6_3)
fitted_param$fit[13] <- "yes"
fitted_param$method[13] <- "nls_lm"
fitted_param$rep[13] <- "3"

#buchanan no lag
b116_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b116_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b116_6_3_index], 
                                         log10nmax = param_estimates$nmax[b116_6_3_index], 
                                         mumax = (param_estimates$mumax[b116_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[14] <-"S12-0116"
fitted_param$pop[14] <-"apc"
fitted_param$model_type[14] <-"buchanan_no_lag"
fitted_param$n0[14] <- summary(b116_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[14] <- summary(b116_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[14] <- summary(b116_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[14] <- "6"
fitted_param$aic[14] <- AIC(b116_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[14] <- BIC(b116_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[14] <- "yes"
fitted_param$method[14] <- "nls_lm"
fitted_param$rep[14] <- "3"

#gompertz
b116_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b116_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b116_6_3_index], 
                                  log10nmax = param_estimates$nmax[b116_6_3_index], 
                                  mumax = (param_estimates$mumax[b116_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b116_6_3_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[15] <-"S12-0116"
fitted_param$pop[15] <-"apc"
fitted_param$model_type[15] <-"gompertz"
fitted_param$n0[15] <- summary(b116_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[15] <- summary(b116_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[15] <- summary(b116_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[15] <- summary(b116_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[15] <- "6"
fitted_param$aic[15] <- AIC(b116_gompertz_nls_lm_6_3)
fitted_param$bic[15] <- BIC(b116_gompertz_nls_lm_6_3)
fitted_param$fit[15] <- "yes"
fitted_param$method[15] <- "nls_lm"
fitted_param$rep[15] <- "3"

#End of S12-0116 model fitting 

##-------------------------S12-0132 Model Fitting, 6C---------------------------

#Rep 1
b132_6_1_index <- which(param_estimates$isolate == "S12-0132" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b132_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b132_6_1, 
                             start=list(
                               log10n0 = param_estimates$n0[b132_6_1_index], 
                               log10nmax = param_estimates$nmax[b132_6_1_index], 
                               mumax = (param_estimates$mumax[b132_6_1_index]*2.303), 
                               lag = param_estimates$lag[b132_6_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[16] <-"S12-0132"
fitted_param$pop[16] <-"apc"
fitted_param$model_type[16] <-"baranyi"
fitted_param$n0[16] <- summary(b132_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[16] <- summary(b132_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[16] <- summary(b132_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[16] <- summary(b132_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[16] <- "6"
fitted_param$aic[16] <- AIC(b132_baranyi_nls_lm_6_1)
fitted_param$bic[16] <- BIC(b132_baranyi_nls_lm_6_1)
fitted_param$fit[16] <- "yes"
fitted_param$method[16] <- "nls_lm"
fitted_param$rep[16] <- "1"

#baranyi no lag
b132_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b132_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b132_6_1_index], 
                                      log10nmax = param_estimates$nmax[b132_6_1_index], 
                                      mumax = (param_estimates$mumax[b132_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[17] <-"S12-0132"
fitted_param$pop[17] <-"apc"
fitted_param$model_type[17] <-"baranyi_no_lag"
fitted_param$n0[17] <- summary(b132_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[17] <- summary(b132_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[17] <- summary(b132_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[17] <- "6"
fitted_param$aic[17] <- AIC(b132_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[17] <- BIC(b132_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[17] <- "yes"
fitted_param$method[17] <- "nls_lm"
fitted_param$rep[17] <- "1"

#buchanan
b132_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b132_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b132_6_1_index], 
                                log10nmax = param_estimates$nmax[b132_6_1_index], 
                                mumax = (param_estimates$mumax[b132_6_1_index]*2.303), 
                                lag = param_estimates$lag[b132_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[18] <-"S12-0132"
fitted_param$pop[18] <-"apc"
fitted_param$model_type[18] <-"buchanan"
fitted_param$n0[18] <- summary(b132_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[18] <- summary(b132_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[18] <- summary(b132_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[18] <- summary(b132_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[18] <- "6"
fitted_param$aic[18] <- AIC(b132_buchanan_nls_lm_6_1)
fitted_param$bic[18] <- BIC(b132_buchanan_nls_lm_6_1)
fitted_param$fit[18] <- "yes"
fitted_param$method[18] <- "nls_lm"
fitted_param$rep[18] <- "1"

#buchanan no lag
b132_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b132_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b132_6_1_index], 
                                       log10nmax = param_estimates$nmax[b132_6_1_index], 
                                       mumax = (param_estimates$mumax[b132_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[19] <-"S12-0132"
fitted_param$pop[19] <-"apc"
fitted_param$model_type[19] <-"buchanan_no_lag"
fitted_param$n0[19] <- summary(b132_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[19] <- summary(b132_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[19] <- summary(b132_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[19] <- "6"
fitted_param$aic[19] <- AIC(b132_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[19] <- BIC(b132_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[19] <- "yes"
fitted_param$method[19] <- "nls_lm"
fitted_param$rep[19] <- "1"

#gompertz
b132_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b132_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b132_6_1_index], 
                                log10nmax = param_estimates$nmax[b132_6_1_index], 
                                mumax = (param_estimates$mumax[b132_6_1_index]*2.303), 
                                lag = param_estimates$lag[b132_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[20] <-"S12-0132"
fitted_param$pop[20] <-"apc"
fitted_param$model_type[20] <-"gompertz"
fitted_param$n0[20] <- summary(b132_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[20] <- summary(b132_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[20] <- summary(b132_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[20] <- summary(b132_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[20] <- "6"
fitted_param$aic[20] <- AIC(b132_gompertz_nls_lm_6_1)
fitted_param$bic[20] <- BIC(b132_gompertz_nls_lm_6_1)
fitted_param$fit[20] <- "yes"
fitted_param$method[20] <- "nls_lm"
fitted_param$rep[20] <- "1"

#Rep 2
b132_6_2_index <- which(param_estimates$isolate == "S12-0132" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b132_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b132_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b132_6_2_index], 
                                 log10nmax = param_estimates$nmax[b132_6_2_index], 
                                 mumax = (param_estimates$mumax[b132_6_2_index]*2.303), 
                                 lag = param_estimates$lag[b132_6_2_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[21] <-"S12-0132"
fitted_param$pop[21] <-"apc"
fitted_param$model_type[21] <-"baranyi"
fitted_param$n0[21] <- summary(b132_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[21] <- summary(b132_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[21] <- summary(b132_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[21] <- summary(b132_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[21] <- "6"
fitted_param$aic[21] <- AIC(b132_baranyi_nls_lm_6_2)
fitted_param$bic[21] <- BIC(b132_baranyi_nls_lm_6_2)
fitted_param$fit[21] <- "yes"
fitted_param$method[21] <- "nls_lm"
fitted_param$rep[21] <- "2"

#baranyi no lag
b132_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b132_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b132_6_2_index], 
                                        log10nmax = param_estimates$nmax[b132_6_2_index], 
                                        mumax = (param_estimates$mumax[b132_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[22] <-"S12-0132"
fitted_param$pop[22] <-"apc"
fitted_param$model_type[22] <-"baranyi_no_lag"
fitted_param$n0[22] <- summary(b132_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[22] <- summary(b132_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[22] <- summary(b132_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[22] <- "6"
fitted_param$aic[22] <- AIC(b132_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[22] <- BIC(b132_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[22] <- "yes"
fitted_param$method[22] <- "nls_lm"
fitted_param$rep[22] <- "2"

#buchanan
b132_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b132_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b132_6_2_index], 
                                  log10nmax = param_estimates$nmax[b132_6_2_index], 
                                  mumax = (param_estimates$mumax[b132_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b132_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[23] <-"S12-0132"
fitted_param$pop[23] <-"apc"
fitted_param$model_type[23] <-"buchanan"
fitted_param$n0[23] <- summary(b132_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[23] <- summary(b132_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[23] <- summary(b132_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[23] <- summary(b132_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[23] <- "6"
fitted_param$aic[23] <- AIC(b132_buchanan_nls_lm_6_2)
fitted_param$bic[23] <- BIC(b132_buchanan_nls_lm_6_2)
fitted_param$fit[23] <- "yes"
fitted_param$method[23] <- "nls_lm"
fitted_param$rep[23] <- "2"


#buchanan no lag
b132_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b132_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b132_6_2_index], 
                                         log10nmax = param_estimates$nmax[b132_6_2_index], 
                                         mumax = (param_estimates$mumax[b132_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[24] <-"S12-0132"
fitted_param$pop[24] <-"apc"
fitted_param$model_type[24] <-"buchanan_no_lag"
fitted_param$n0[24] <- summary(b132_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[24] <- summary(b132_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[24] <- summary(b132_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[24] <- "6"
fitted_param$aic[24] <- AIC(b132_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[24] <- BIC(b132_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[24] <- "yes"
fitted_param$method[24] <- "nls_lm"
fitted_param$rep[24] <- "2"

#gompertz
b132_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b132_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b132_6_2_index], 
                                  log10nmax = param_estimates$nmax[b132_6_2_index], 
                                  mumax = (param_estimates$mumax[b132_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b132_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[25] <-"S12-0132"
fitted_param$pop[25] <-"apc"
fitted_param$model_type[25] <-"gompertz"
fitted_param$n0[25] <- summary(b132_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[25] <- summary(b132_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[25] <- summary(b132_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[25] <- summary(b132_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[25] <- "6"
fitted_param$aic[25] <- AIC(b132_gompertz_nls_lm_6_2)
fitted_param$bic[25] <- BIC(b132_gompertz_nls_lm_6_2)
fitted_param$fit[25] <- "yes"
fitted_param$method[25] <- "nls_lm"
fitted_param$rep[25] <- "2"


#Rep 3
b132_6_3_index <- which(param_estimates$isolate == "S12-0132" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b132_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b132_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b132_6_3_index], 
                                 log10nmax = param_estimates$nmax[b132_6_3_index], 
                                 mumax = (param_estimates$mumax[b132_6_3_index]*2.303), 
                                 lag = param_estimates$lag[b132_6_3_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[26] <-"S12-0132"
fitted_param$pop[26] <-"apc"
fitted_param$model_type[26] <-"baranyi"
fitted_param$n0[26] <- summary(b132_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[26] <- summary(b132_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[26] <- summary(b132_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[26] <- summary(b132_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[26] <- "6"
fitted_param$aic[26] <- AIC(b132_baranyi_nls_lm_6_3)
fitted_param$bic[26] <- BIC(b132_baranyi_nls_lm_6_3)
fitted_param$fit[26] <- "yes"
fitted_param$method[26] <- "nls_lm"
fitted_param$rep[26] <- "3"

#baranyi no lag
b132_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b132_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b132_6_3_index], 
                                        log10nmax = param_estimates$nmax[b132_6_3_index], 
                                        mumax = (param_estimates$mumax[b132_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[27] <-"S12-0132"
fitted_param$pop[27] <-"apc"
fitted_param$model_type[27] <-"baranyi_no_lag"
fitted_param$n0[27] <- summary(b132_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[27] <- summary(b132_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[27] <- summary(b132_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[27] <- "6"
fitted_param$aic[27] <- AIC(b132_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[27] <- BIC(b132_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[27] <- "yes"
fitted_param$method[27] <- "nls_lm"
fitted_param$rep[27] <- "3"


#buchanan
b132_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b132_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b132_6_3_index], 
                                  log10nmax = param_estimates$nmax[b132_6_3_index], 
                                  mumax = (param_estimates$mumax[b132_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b132_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[28] <-"S12-0132"
fitted_param$pop[28] <-"apc"
fitted_param$model_type[28] <-"buchanan"
fitted_param$n0[28] <- summary(b132_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[28] <- summary(b132_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[28] <- summary(b132_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[28] <- summary(b132_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[28] <- "6"
fitted_param$aic[28] <- AIC(b132_buchanan_nls_lm_6_3)
fitted_param$bic[28] <- BIC(b132_buchanan_nls_lm_6_3)
fitted_param$fit[28] <- "no"
fitted_param$method[28] <- "nls_lm"
fitted_param$rep[28] <- "3"

#buchanan no lag
b132_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b132_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b132_6_3_index], 
                                         log10nmax = param_estimates$nmax[b132_6_3_index], 
                                         mumax = (param_estimates$mumax[b132_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[29] <-"S12-0132"
fitted_param$pop[29] <-"apc"
fitted_param$model_type[29] <-"buchanan_no_lag"
fitted_param$n0[29] <- summary(b132_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[29] <- summary(b132_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[29] <- summary(b132_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[29] <- "6"
fitted_param$aic[29] <- AIC(b132_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[29] <- BIC(b132_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[29] <- "yes"
fitted_param$method[29] <- "nls_lm"
fitted_param$rep[29] <- "3"

#gompertz
b132_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b132_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b132_6_3_index], 
                                  log10nmax = param_estimates$nmax[b132_6_3_index], 
                                  mumax = (param_estimates$mumax[b132_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b132_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[30] <-"S12-0132"
fitted_param$pop[30] <-"apc"
fitted_param$model_type[30] <-"gompertz"
fitted_param$n0[30] <- summary(b132_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[30] <- summary(b132_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[30] <- summary(b132_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[30] <- summary(b132_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[30] <- "6"
fitted_param$aic[30] <- AIC(b132_gompertz_nls_lm_6_3)
fitted_param$bic[30] <- BIC(b132_gompertz_nls_lm_6_3)
fitted_param$fit[30] <- "yes"
fitted_param$method[30] <- "nls_lm"
fitted_param$rep[30] <- "3"


#End of S12-0132 model fitting 

##------------------------S12-0141 Model Fitting, 6C----------------------------

#Rep 1
b141_6_1_index <- which(param_estimates$isolate == "S12-0141" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b141_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b141_6_1, 
                             start=list(
                               log10n0 = param_estimates$n0[b141_6_1_index], 
                               log10nmax = param_estimates$nmax[b141_6_1_index], 
                               mumax = (param_estimates$mumax[b141_6_1_index]*2.303), 
                               lag = param_estimates$lag[b141_6_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[31] <-"S12-0141"
fitted_param$pop[31] <-"apc"
fitted_param$model_type[31] <-"baranyi"
fitted_param$n0[31] <- summary(b141_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[31] <- summary(b141_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[31] <- summary(b141_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[31] <- summary(b141_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[31] <- "6"
fitted_param$aic[31] <- AIC(b141_baranyi_nls_lm_6_1)
fitted_param$bic[31] <- BIC(b141_baranyi_nls_lm_6_1)
fitted_param$fit[31] <- "yes"
fitted_param$method[31] <- "nls_lm"
fitted_param$rep[31] <- "1"


#baranyi no lag
b141_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b141_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b141_6_1_index], 
                                      log10nmax = param_estimates$nmax[b141_6_1_index], 
                                      mumax = (param_estimates$mumax[b141_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[32] <-"S12-0141"
fitted_param$pop[32] <-"apc"
fitted_param$model_type[32] <-"baranyi_no_lag"
fitted_param$n0[32] <- summary(b141_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[32] <- summary(b141_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[32] <- summary(b141_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[32] <- "6"
fitted_param$aic[32] <- AIC(b141_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[32] <- BIC(b141_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[32] <- "yes"
fitted_param$method[32] <- "nls_lm"
fitted_param$rep[32] <- "1"

#buchanan
b141_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b141_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b141_6_1_index], 
                                log10nmax = param_estimates$nmax[b141_6_1_index], 
                                mumax = (param_estimates$mumax[b141_6_1_index]*2.303), 
                                lag = param_estimates$lag[b141_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[33] <-"S12-0141"
fitted_param$pop[33] <-"apc"
fitted_param$model_type[33] <-"buchanan"
fitted_param$n0[33] <- summary(b141_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[33] <- summary(b141_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[33] <- summary(b141_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[33] <- summary(b141_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[33] <- "6"
fitted_param$aic[33] <- AIC(b141_buchanan_nls_lm_6_1)
fitted_param$bic[33] <- BIC(b141_buchanan_nls_lm_6_1)
fitted_param$fit[33] <- "yes"
fitted_param$method[33] <- "nls_lm"
fitted_param$rep[33] <- "1"


#buchanan no lag
b141_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b141_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b141_6_1_index], 
                                       log10nmax = param_estimates$nmax[b141_6_1_index], 
                                       mumax = (param_estimates$mumax[b141_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[34] <-"S12-0141"
fitted_param$pop[34] <-"apc"
fitted_param$model_type[34] <-"buchanan_no_lag"
fitted_param$n0[34] <- summary(b141_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[34] <- summary(b141_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[34] <- summary(b141_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[34] <- "6"
fitted_param$aic[34] <- AIC(b141_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[34] <- BIC(b141_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[34] <- "yes"
fitted_param$method[34] <- "nls_lm"
fitted_param$rep[34] <- "1"


#gompertz
b141_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b141_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b141_6_1_index], 
                                log10nmax = param_estimates$nmax[b141_6_1_index], 
                                mumax = (param_estimates$mumax[b141_6_1_index]*2.303), 
                                lag = param_estimates$lag[b141_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[35] <-"S12-0141"
fitted_param$pop[35] <-"apc"
fitted_param$model_type[35] <-"gompertz"
fitted_param$n0[35] <- summary(b141_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[35] <- summary(b141_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[35] <- summary(b141_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[35] <- summary(b141_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[35] <- "6"
fitted_param$aic[35] <- AIC(b141_gompertz_nls_lm_6_1)
fitted_param$bic[35] <- BIC(b141_gompertz_nls_lm_6_1)
fitted_param$fit[35] <- "yes"
fitted_param$method[35] <- "nls_lm"
fitted_param$rep[35] <- "1"



#Rep 2
b141_6_2_index <- which(param_estimates$isolate == "S12-0141" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b141_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b141_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b141_6_2_index], 
                                 log10nmax = param_estimates$nmax[b141_6_2_index], 
                                 mumax = (param_estimates$mumax[b141_6_2_index]*2.303), 
                                 lag = param_estimates$lag[b141_6_2_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[36] <-"S12-0141"
fitted_param$pop[36] <-"apc"
fitted_param$model_type[36] <-"baranyi"
fitted_param$n0[36] <- summary(b141_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[36] <- summary(b141_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[36] <- summary(b141_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[36] <- summary(b141_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[36] <- "6"
fitted_param$aic[36] <- AIC(b141_baranyi_nls_lm_6_2)
fitted_param$bic[36] <- BIC(b141_baranyi_nls_lm_6_2)
fitted_param$fit[36] <- "yes"
fitted_param$method[36] <- "nls_lm"
fitted_param$rep[36] <- "2"


#baranyi no lag
b141_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b141_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b141_6_2_index], 
                                        log10nmax = param_estimates$nmax[b141_6_2_index], 
                                        mumax = (param_estimates$mumax[b141_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[37] <-"S12-0141"
fitted_param$pop[37] <-"apc"
fitted_param$model_type[37] <-"baranyi_no_lag"
fitted_param$n0[37] <- summary(b141_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[37] <- summary(b141_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[37] <- summary(b141_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[37] <- "6"
fitted_param$aic[37] <- AIC(b141_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[37] <- BIC(b141_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[37] <- "yes"
fitted_param$method[37] <- "nls_lm"
fitted_param$rep[37] <- "2"

#buchanan
b141_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b141_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b141_6_2_index], 
                                  log10nmax = param_estimates$nmax[b141_6_2_index], 
                                  mumax = (param_estimates$mumax[b141_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b141_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[38] <-"S12-0141"
fitted_param$pop[38] <-"apc"
fitted_param$model_type[38] <-"buchanan"
fitted_param$n0[38] <- summary(b141_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[38] <- summary(b141_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[38] <- summary(b141_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[38] <- summary(b141_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[38] <- "6"
fitted_param$aic[38] <- AIC(b141_buchanan_nls_lm_6_2)
fitted_param$bic[38] <- BIC(b141_buchanan_nls_lm_6_2)
fitted_param$fit[38] <- "yes"
fitted_param$method[38] <- "nls_lm"
fitted_param$rep[38] <- "2"

#buchanan no lag
b141_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b141_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b141_6_2_index], 
                                         log10nmax = param_estimates$nmax[b141_6_2_index], 
                                         mumax = (param_estimates$mumax[b141_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[39] <-"S12-0141"
fitted_param$pop[39] <-"apc"
fitted_param$model_type[39] <-"buchanan_no_lag"
fitted_param$n0[39] <- summary(b141_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[39] <- summary(b141_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[39] <- summary(b141_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[39] <- "6"
fitted_param$aic[39] <- AIC(b141_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[39] <- BIC(b141_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[39] <- "yes"
fitted_param$method[39] <- "nls_lm"
fitted_param$rep[39] <- "2"

#gompertz
b141_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b141_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b141_6_2_index], 
                                  log10nmax = param_estimates$nmax[b141_6_2_index], 
                                  mumax = (param_estimates$mumax[b141_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b141_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[40] <-"S12-0141"
fitted_param$pop[40] <-"apc"
fitted_param$model_type[40] <-"gompertz"
fitted_param$n0[40] <- summary(b141_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[40] <- summary(b141_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[40] <- summary(b141_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[40] <- summary(b141_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[40] <- "6"
fitted_param$aic[40] <- AIC(b141_gompertz_nls_lm_6_2)
fitted_param$bic[40] <- BIC(b141_gompertz_nls_lm_6_2)
fitted_param$fit[40] <- "yes"
fitted_param$method[40] <- "nls_lm"
fitted_param$rep[40] <- "2"

#Rep 3
b141_6_3_index <- which(param_estimates$isolate == "S12-0141" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b141_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b141_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b141_6_3_index], 
                                 log10nmax = param_estimates$nmax[b141_6_3_index], 
                                 mumax = (param_estimates$mumax[b141_6_3_index]*2.303), 
                                 lag = param_estimates$lag[b141_6_3_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[41] <-"S12-0141"
fitted_param$pop[41] <-"apc"
fitted_param$model_type[41] <-"baranyi"
fitted_param$n0[41] <- summary(b141_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[41] <- summary(b141_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[41] <- summary(b141_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[41] <- summary(b141_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[41] <- "6"
fitted_param$aic[41] <- AIC(b141_baranyi_nls_lm_6_3)
fitted_param$bic[41] <- BIC(b141_baranyi_nls_lm_6_3)
fitted_param$fit[41] <- "yes"
fitted_param$method[41] <- "nls_lm"
fitted_param$rep[41] <- "3"

#baranyi no lag
b141_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b141_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b141_6_3_index], 
                                        log10nmax = param_estimates$nmax[b141_6_3_index], 
                                        mumax = (param_estimates$mumax[b141_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[42] <-"S12-0141"
fitted_param$pop[42] <-"apc"
fitted_param$model_type[42] <-"baranyi_no_lag"
fitted_param$n0[42] <- summary(b141_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[42] <- summary(b141_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[42] <- summary(b141_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[42] <- "6"
fitted_param$aic[42] <- AIC(b141_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[42] <- BIC(b141_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[42] <- "yes"
fitted_param$method[42] <- "nls_lm"
fitted_param$rep[42] <- "3"

#buchanan
b141_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b141_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b141_6_3_index], 
                                  log10nmax = param_estimates$nmax[b141_6_3_index], 
                                  mumax = (param_estimates$mumax[b141_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b141_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[43] <-"S12-0141"
fitted_param$pop[43] <-"apc"
fitted_param$model_type[43] <-"buchanan"
fitted_param$n0[43] <- summary(b141_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[43] <- summary(b141_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[43] <- summary(b141_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[43] <- summary(b141_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[43] <- "6"
fitted_param$aic[43] <- AIC(b141_buchanan_nls_lm_6_3)
fitted_param$bic[43] <- BIC(b141_buchanan_nls_lm_6_3)
fitted_param$fit[43] <- "yes"
fitted_param$method[43] <- "nls_lm"
fitted_param$rep[43] <- "3"

#buchanan no lag
b141_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b141_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b141_6_3_index], 
                                         log10nmax = param_estimates$nmax[b141_6_3_index], 
                                         mumax = (param_estimates$mumax[b141_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[44] <-"S12-0141"
fitted_param$pop[44] <-"apc"
fitted_param$model_type[44] <-"buchanan_no_lag"
fitted_param$n0[44] <- summary(b141_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[44] <- summary(b141_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[44] <- summary(b141_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[44] <- "6"
fitted_param$aic[44] <- AIC(b141_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[44] <- BIC(b141_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[44] <- "yes"
fitted_param$method[44] <- "nls_lm"
fitted_param$rep[44] <- "3"

#gompertz
b141_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b141_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b141_6_3_index], 
                                  log10nmax = param_estimates$nmax[b141_6_3_index], 
                                  mumax = (param_estimates$mumax[b141_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b141_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[45] <-"S12-0141"
fitted_param$pop[45] <-"apc"
fitted_param$model_type[45] <-"gompertz"
fitted_param$n0[45] <- summary(b141_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[45] <- summary(b141_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[45] <- summary(b141_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[45] <- summary(b141_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[45] <- "6"
fitted_param$aic[45] <- AIC(b141_gompertz_nls_lm_6_3)
fitted_param$bic[45] <- BIC(b141_gompertz_nls_lm_6_3)
fitted_param$fit[45] <- "yes"
fitted_param$method[45] <- "nls_lm"
fitted_param$rep[45] <- "3"

#End of S12-0141 model fitting 

##------------------------S12-0166 Model Fitting, 6C----------------------------

#Rep 1
b166_6_1_index <- which(param_estimates$isolate == "S12-0166" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b166_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b166_6_1, 
                             start=list(
                               log10n0 = param_estimates$n0[b166_6_1_index], 
                               log10nmax = param_estimates$nmax[b166_6_1_index], 
                               mumax = (param_estimates$mumax[b166_6_1_index]*2.303), 
                               lag = param_estimates$lag[b166_6_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[46] <-"S12-0166"
fitted_param$pop[46] <-"apc"
fitted_param$model_type[46] <-"baranyi"
fitted_param$n0[46] <- summary(b166_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[46] <- summary(b166_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[46] <- summary(b166_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[46] <- summary(b166_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[46] <- "6"
fitted_param$aic[46] <- AIC(b166_baranyi_nls_lm_6_1)
fitted_param$bic[46] <- BIC(b166_baranyi_nls_lm_6_1)
fitted_param$fit[46] <- "yes"
fitted_param$method[46] <- "nls_lm"
fitted_param$rep[46] <- "1"

#baranyi no lag
b166_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b166_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b166_6_1_index], 
                                      log10nmax = param_estimates$nmax[b166_6_1_index], 
                                      mumax = (param_estimates$mumax[b166_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[47] <-"S12-0166"
fitted_param$pop[47] <-"apc"
fitted_param$model_type[47] <-"baranyi_no_lag"
fitted_param$n0[47] <- summary(b166_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[47] <- summary(b166_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[47] <- summary(b166_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[47] <- "6"
fitted_param$aic[47] <- AIC(b166_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[47] <- BIC(b166_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[47] <- "yes"
fitted_param$method[47] <- "nls_lm"
fitted_param$rep[47] <- "1"

#buchanan
b166_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b166_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b166_6_1_index], 
                                log10nmax = param_estimates$nmax[b166_6_1_index], 
                                mumax = (param_estimates$mumax[b166_6_1_index]*2.303), 
                                lag = param_estimates$lag[b166_6_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[48] <-"S12-0166"
fitted_param$pop[48] <-"apc"
fitted_param$model_type[48] <-"buchanan"
fitted_param$n0[48] <- summary(b166_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[48] <- summary(b166_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[48] <- summary(b166_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[48] <- summary(b166_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[48] <- "6"
fitted_param$aic[48] <- AIC(b166_buchanan_nls_lm_6_1)
fitted_param$bic[48] <- BIC(b166_buchanan_nls_lm_6_1)
fitted_param$fit[48] <- "yes"
fitted_param$method[48] <- "nls_lm"
fitted_param$rep[48] <- "1"

#buchanan no lag
b166_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b166_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b166_6_1_index], 
                                       log10nmax = param_estimates$nmax[b166_6_1_index], 
                                       mumax = (param_estimates$mumax[b166_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[49] <-"S12-0166"
fitted_param$pop[49] <-"apc"
fitted_param$model_type[49] <-"buchanan_no_lag"
fitted_param$n0[49] <- summary(b166_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[49] <- summary(b166_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[49] <- summary(b166_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[49] <- "6"
fitted_param$aic[49] <- AIC(b166_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[49] <- BIC(b166_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[49] <- "yes"
fitted_param$method[49] <- "nls_lm"
fitted_param$rep[49] <- "1"

#gompertz
b166_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b166_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b166_6_1_index], 
                                log10nmax = param_estimates$nmax[b166_6_1_index], 
                                mumax = (param_estimates$mumax[b166_6_1_index]*2.303),
                                lag = param_estimates$lag[b166_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[50] <-"S12-0166"
fitted_param$pop[50] <-"apc"
fitted_param$model_type[50] <-"gompertz"
fitted_param$n0[50] <- summary(b166_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[50] <- summary(b166_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[50] <- summary(b166_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[50] <- summary(b166_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[50] <- "6"
fitted_param$aic[50] <- AIC(b166_gompertz_nls_lm_6_1)
fitted_param$bic[50] <- BIC(b166_gompertz_nls_lm_6_1)
fitted_param$fit[50] <- "yes"
fitted_param$method[50] <- "nls_lm"
fitted_param$rep[50] <- "1"

#Rep 2
b166_6_2_index <- which(param_estimates$isolate == "S12-0166" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b166_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b166_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b166_6_2_index], 
                                 log10nmax = param_estimates$nmax[b166_6_2_index], 
                                 mumax = (param_estimates$mumax[b166_6_2_index]*2.303), 
                                 lag = param_estimates$lag[b166_6_2_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[51] <-"S12-0166"
fitted_param$pop[51] <-"apc"
fitted_param$model_type[51] <-"baranyi"
fitted_param$n0[51] <- summary(b166_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[51] <- summary(b166_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[51] <- summary(b166_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[51] <- summary(b166_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[51] <- "6"
fitted_param$aic[51] <- AIC(b166_baranyi_nls_lm_6_2)
fitted_param$bic[51] <- BIC(b166_baranyi_nls_lm_6_2)
fitted_param$fit[51] <- "no"
fitted_param$method[51] <- "nls_lm"
fitted_param$rep[51] <- "2"

#baranyi no lag
b166_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b166_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b166_6_2_index], 
                                        log10nmax = param_estimates$nmax[b166_6_2_index], 
                                        mumax = (param_estimates$mumax[b166_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[52] <-"S12-0166"
fitted_param$pop[52] <-"apc"
fitted_param$model_type[52] <-"baranyi_no_lag"
fitted_param$n0[52] <- summary(b166_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[52] <- summary(b166_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[52] <- summary(b166_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[52] <- "6"
fitted_param$aic[52] <- AIC(b166_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[52] <- BIC(b166_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[52] <- "no"
fitted_param$method[52] <- "nls_lm"
fitted_param$rep[52] <- "2"

#buchanan
b166_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b166_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b166_6_2_index], 
                                  log10nmax = param_estimates$nmax[b166_6_2_index], 
                                  mumax = (param_estimates$mumax[b166_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b166_6_2_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[53] <-"S12-0166"
fitted_param$pop[53] <-"apc"
fitted_param$model_type[53] <-"buchanan"
fitted_param$n0[53] <- summary(b166_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[53] <- summary(b166_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[53] <- summary(b166_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[53] <- summary(b166_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[53] <- "6"
fitted_param$aic[53] <- AIC(b166_buchanan_nls_lm_6_2)
fitted_param$bic[53] <- BIC(b166_buchanan_nls_lm_6_2)
fitted_param$fit[53] <- "yes"
fitted_param$method[53] <- "nls_lm"
fitted_param$rep[53] <- "2"

#buchanan no lag
b166_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b166_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b166_6_2_index], 
                                         log10nmax = param_estimates$nmax[b166_6_2_index], 
                                         mumax = (param_estimates$mumax[b166_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[54] <-"S12-0166"
fitted_param$pop[54] <-"apc"
fitted_param$model_type[54] <-"buchanan_no_lag"
fitted_param$n0[54] <- summary(b166_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[54] <- summary(b166_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[54] <- summary(b166_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[54] <- "6"
fitted_param$aic[54] <- AIC(b166_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[54] <- BIC(b166_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[54] <- "yes"
fitted_param$method[54] <- "nls_lm"
fitted_param$rep[54] <- "2"

#gompertz
b166_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b166_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b166_6_2_index], 
                                  log10nmax = param_estimates$nmax[b166_6_2_index], 
                                  mumax = (param_estimates$mumax[b166_6_2_index]*2.303),
                                  lag = param_estimates$lag[b166_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[55] <-"S12-0166"
fitted_param$pop[55] <-"apc"
fitted_param$model_type[55] <-"gompertz"
fitted_param$n0[55] <- summary(b166_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[55] <- summary(b166_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[55] <- summary(b166_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[55] <- summary(b166_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[55] <- "6"
fitted_param$aic[55] <- AIC(b166_gompertz_nls_lm_6_2)
fitted_param$bic[55] <- BIC(b166_gompertz_nls_lm_6_2)
fitted_param$fit[55] <- "yes"
fitted_param$method[55] <- "nls_lm"
fitted_param$rep[55] <- "2"

#Rep 3
b166_6_3_index <- which(param_estimates$isolate == "S12-0166" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b166_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b166_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b166_6_3_index], 
                                 log10nmax = param_estimates$nmax[b166_6_3_index], 
                                 mumax = (param_estimates$mumax[b166_6_3_index]*2.303), 
                                 lag = param_estimates$lag[b166_6_3_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[56] <-"S12-0166"
fitted_param$pop[56] <-"apc"
fitted_param$model_type[56] <-"baranyi"
fitted_param$n0[56] <- summary(b166_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[56] <- summary(b166_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[56] <- summary(b166_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[56] <- summary(b166_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[56] <- "6"
fitted_param$aic[56] <- AIC(b166_baranyi_nls_lm_6_3)
fitted_param$bic[56] <- BIC(b166_baranyi_nls_lm_6_3)
fitted_param$fit[56] <- "yes"
fitted_param$method[56] <- "nls_lm"
fitted_param$rep[56] <- "3"

#baranyi no lag
b166_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b166_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b166_6_3_index], 
                                        log10nmax = param_estimates$nmax[b166_6_3_index], 
                                        mumax = (param_estimates$mumax[b166_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[57] <-"S12-0166"
fitted_param$pop[57] <-"apc"
fitted_param$model_type[57] <-"baranyi_no_lag"
fitted_param$n0[57] <- summary(b166_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[57] <- summary(b166_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[57] <- summary(b166_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[57] <- "6"
fitted_param$aic[57] <- AIC(b166_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[57] <- BIC(b166_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[57] <- "yes"
fitted_param$method[57] <- "nls_lm"
fitted_param$rep[57] <- "3"

#buchanan
b166_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b166_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b166_6_3_index], 
                                  log10nmax = param_estimates$nmax[b166_6_3_index], 
                                  mumax = (param_estimates$mumax[b166_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b166_6_3_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[58] <-"S12-0166"
fitted_param$pop[58] <-"apc"
fitted_param$model_type[58] <-"buchanan"
fitted_param$n0[58] <- summary(b166_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[58] <- summary(b166_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[58] <- summary(b166_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[58] <- summary(b166_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[58] <- "6"
fitted_param$aic[58] <- AIC(b166_buchanan_nls_lm_6_3)
fitted_param$bic[58] <- BIC(b166_buchanan_nls_lm_6_3)
fitted_param$fit[58] <- "yes"
fitted_param$method[58] <- "nls_lm"
fitted_param$rep[58] <- "3"

#buchanan no lag
b166_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b166_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b166_6_3_index], 
                                         log10nmax = param_estimates$nmax[b166_6_3_index], 
                                         mumax = (param_estimates$mumax[b166_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[59] <-"S12-0166"
fitted_param$pop[59] <-"apc"
fitted_param$model_type[59] <-"buchanan_no_lag"
fitted_param$n0[59] <- summary(b166_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[59] <- summary(b166_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[59] <- summary(b166_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[59] <- "6"
fitted_param$aic[59] <- AIC(b166_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[59] <- BIC(b166_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[59] <- "yes"
fitted_param$method[59] <- "nls_lm"
fitted_param$rep[59] <- "3"

#gompertz
b166_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b166_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b166_6_3_index], 
                                  log10nmax = param_estimates$nmax[b166_6_3_index], 
                                  mumax = (param_estimates$mumax[b166_6_3_index]*2.303),
                                  lag = param_estimates$lag[b166_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[60] <-"S12-0166"
fitted_param$pop[60] <-"apc"
fitted_param$model_type[60] <-"gompertz"
fitted_param$n0[60] <- summary(b166_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[60] <- summary(b166_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[60] <- summary(b166_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[60] <- summary(b166_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[60] <- "6"
fitted_param$aic[60] <- AIC(b166_gompertz_nls_lm_6_3)
fitted_param$bic[60] <- BIC(b166_gompertz_nls_lm_6_3)
fitted_param$fit[60] <- "yes"
fitted_param$method[60] <- "nls_lm"
fitted_param$rep[60] <- "3"

#End of S12-0166 model fitting 

##------------------------S12-0180 Model Fitting, 6C----------------------------

#Rep 1
b180_6_1_index <- which(param_estimates$isolate == "S12-0180" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b180_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b180_6_1, 
                             start=list(
                               log10n0 = param_estimates$n0[b180_6_1_index], 
                               log10nmax = param_estimates$nmax[b180_6_1_index], 
                               mumax = (param_estimates$mumax[b180_6_1_index]*2.303),
                               lag = param_estimates$lag[b180_6_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[61] <-"S12-0180"
fitted_param$pop[61] <-"apc"
fitted_param$model_type[61] <-"baranyi"
fitted_param$n0[61] <- summary(b180_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[61] <- summary(b180_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[61] <- summary(b180_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[61] <- summary(b180_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[61] <- "6"
fitted_param$aic[61] <- AIC(b180_baranyi_nls_lm_6_1)
fitted_param$bic[61] <- BIC(b180_baranyi_nls_lm_6_1)
fitted_param$fit[61] <- "yes"
fitted_param$method[61] <- "nls_lm"
fitted_param$rep[61] <- "1"

#baranyi no lag
b180_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b180_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b180_6_1_index], 
                                      log10nmax = param_estimates$nmax[b180_6_1_index], 
                                      mumax = (param_estimates$mumax[b180_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[62] <-"S12-0180"
fitted_param$pop[62] <-"apc"
fitted_param$model_type[62] <-"baranyi_no_lag"
fitted_param$n0[62] <- summary(b180_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[62] <- summary(b180_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[62] <- summary(b180_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[62] <- "6"
fitted_param$aic[62] <- AIC(b180_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[62] <- BIC(b180_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[62] <- "yes"
fitted_param$method[62] <- "nls_lm"
fitted_param$rep[62] <- "1"

#buchanan
b180_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b180_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b180_6_1_index], 
                                log10nmax = param_estimates$nmax[b180_6_1_index], 
                                mumax = (param_estimates$mumax[b180_6_1_index]*2.303),
                                lag = param_estimates$lag[b180_6_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[63] <-"S12-0180"
fitted_param$pop[63] <-"apc"
fitted_param$model_type[63] <-"buchanan"
fitted_param$n0[63] <- summary(b180_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[63] <- summary(b180_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[63] <- summary(b180_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[63] <- summary(b180_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[63] <- "6"
fitted_param$aic[63] <- AIC(b180_buchanan_nls_lm_6_1)
fitted_param$bic[63] <- BIC(b180_buchanan_nls_lm_6_1)
fitted_param$fit[63] <- "yes"
fitted_param$method[63] <- "nls_lm"
fitted_param$rep[63] <- "1"

#buchanan no lag
b180_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b180_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b180_6_1_index], 
                                       log10nmax = param_estimates$nmax[b180_6_1_index], 
                                       mumax = (param_estimates$mumax[b180_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[64] <-"S12-0180"
fitted_param$pop[64] <-"apc"
fitted_param$model_type[64] <-"buchanan_no_lag"
fitted_param$n0[64] <- summary(b180_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[64] <- summary(b180_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[64] <- summary(b180_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[64] <- "6"
fitted_param$aic[64] <- AIC(b180_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[64] <- BIC(b180_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[64] <- "yes"
fitted_param$method[64] <- "nls_lm"
fitted_param$rep[64] <- "1"

#gompertz
b180_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b180_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b180_6_1_index], 
                                log10nmax = param_estimates$nmax[b180_6_1_index], 
                                mumax = (param_estimates$mumax[b180_6_1_index]*2.303),
                                lag = param_estimates$lag[b180_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[65] <-"S12-0180"
fitted_param$pop[65] <-"apc"
fitted_param$model_type[65] <-"gompertz"
fitted_param$n0[65] <- summary(b180_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[65] <- summary(b180_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[65] <- summary(b180_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[65] <- summary(b180_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[65] <- "6"
fitted_param$aic[65] <- AIC(b180_gompertz_nls_lm_6_1)
fitted_param$bic[65] <- BIC(b180_gompertz_nls_lm_6_1)
fitted_param$fit[65] <- "yes"
fitted_param$method[65] <- "nls_lm"
fitted_param$rep[65] <- "1"

#Rep 2
b180_6_2_index <- which(param_estimates$isolate == "S12-0180" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b180_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b180_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b180_6_2_index], 
                                 log10nmax = param_estimates$nmax[b180_6_2_index], 
                                 mumax = (param_estimates$mumax[b180_6_2_index]*2.303),
                                 lag = param_estimates$lag[b180_6_2_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[66] <-"S12-0180"
fitted_param$pop[66] <-"apc"
fitted_param$model_type[66] <-"baranyi"
fitted_param$n0[66] <- summary(b180_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[66] <- summary(b180_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[66] <- summary(b180_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[66] <- summary(b180_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[66] <- "6"
fitted_param$aic[66] <- AIC(b180_baranyi_nls_lm_6_2)
fitted_param$bic[66] <- BIC(b180_baranyi_nls_lm_6_2)
fitted_param$fit[66] <- "no"
fitted_param$method[66] <- "nls_lm"
fitted_param$rep[66] <- "2"

#baranyi no lag
b180_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b180_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b180_6_2_index], 
                                        log10nmax = param_estimates$nmax[b180_6_2_index], 
                                        mumax = (param_estimates$mumax[b180_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[67] <-"S12-0180"
fitted_param$pop[67] <-"apc"
fitted_param$model_type[67] <-"baranyi_no_lag"
fitted_param$n0[67] <- summary(b180_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[67] <- summary(b180_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[67] <- summary(b180_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[67] <- "6"
fitted_param$aic[67] <- AIC(b180_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[67] <- BIC(b180_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[67] <- "yes"
fitted_param$method[67] <- "nls_lm"
fitted_param$rep[67] <- "2"

#buchanan
b180_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b180_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b180_6_2_index], 
                                  log10nmax = param_estimates$nmax[b180_6_2_index], 
                                  mumax = (param_estimates$mumax[b180_6_2_index]*2.303),
                                  lag = param_estimates$lag[b180_6_2_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[68] <-"S12-0180"
fitted_param$pop[68] <-"apc"
fitted_param$model_type[68] <-"buchanan"
fitted_param$n0[68] <- summary(b180_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[68] <- summary(b180_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[68] <- summary(b180_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[68] <- summary(b180_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[68] <- "6"
fitted_param$aic[68] <- AIC(b180_buchanan_nls_lm_6_2)
fitted_param$bic[68] <- BIC(b180_buchanan_nls_lm_6_2)
fitted_param$fit[68] <- "yes"
fitted_param$method[68] <- "nls_lm"
fitted_param$rep[68] <- "2"

#buchanan no lag
b180_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b180_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b180_6_2_index], 
                                         log10nmax = param_estimates$nmax[b180_6_2_index], 
                                         mumax = (param_estimates$mumax[b180_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[69] <-"S12-0180"
fitted_param$pop[69] <-"apc"
fitted_param$model_type[69] <-"buchanan_no_lag"
fitted_param$n0[69] <- summary(b180_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[69] <- summary(b180_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[69] <- summary(b180_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[69] <- "6"
fitted_param$aic[69] <- AIC(b180_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[69] <- BIC(b180_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[69] <- "yes"
fitted_param$method[69] <- "nls_lm"
fitted_param$rep[69] <- "2"

#gompertz
b180_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b180_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b180_6_2_index], 
                                  log10nmax = param_estimates$nmax[b180_6_2_index], 
                                  mumax = (param_estimates$mumax[b180_6_2_index]*2.303),
                                  lag = param_estimates$lag[b180_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[70] <-"S12-0180"
fitted_param$pop[70] <-"apc"
fitted_param$model_type[70] <-"gompertz"
fitted_param$n0[70] <- summary(b180_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[70] <- summary(b180_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[70] <- summary(b180_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[70] <- summary(b180_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[70] <- "6"
fitted_param$aic[70] <- AIC(b180_gompertz_nls_lm_6_2)
fitted_param$bic[70] <- BIC(b180_gompertz_nls_lm_6_2)
fitted_param$fit[70] <- "yes"
fitted_param$method[70] <- "nls_lm"
fitted_param$rep[70] <- "2"

#Rep 3
b180_6_3_index <- which(param_estimates$isolate == "S12-0180" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b180_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b180_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b180_6_3_index], 
                                 log10nmax = param_estimates$nmax[b180_6_3_index], 
                                 mumax = (param_estimates$mumax[b180_6_3_index]*2.303),
                                 lag = param_estimates$lag[b180_6_3_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[71] <-"S12-0180"
fitted_param$pop[71] <-"apc"
fitted_param$model_type[71] <-"baranyi"
fitted_param$n0[71] <- summary(b180_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[71] <- summary(b180_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[71] <- summary(b180_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[71] <- summary(b180_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[71] <- "6"
fitted_param$aic[71] <- AIC(b180_baranyi_nls_lm_6_3)
fitted_param$bic[71] <- BIC(b180_baranyi_nls_lm_6_3)
fitted_param$fit[71] <- "yes"
fitted_param$method[71] <- "nls_lm"
fitted_param$rep[71] <- "3"

#baranyi no lag
b180_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b180_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b180_6_3_index], 
                                        log10nmax = param_estimates$nmax[b180_6_3_index], 
                                        mumax = (param_estimates$mumax[b180_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[72] <-"S12-0180"
fitted_param$pop[72] <-"apc"
fitted_param$model_type[72] <-"baranyi_no_lag"
fitted_param$n0[72] <- summary(b180_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[72] <- summary(b180_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[72] <- summary(b180_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[72] <- "6"
fitted_param$aic[72] <- AIC(b180_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[72] <- BIC(b180_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[72] <- "no"
fitted_param$method[72] <- "nls_lm"
fitted_param$rep[72] <- "3"

#buchanan
b180_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b180_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b180_6_3_index], 
                                  log10nmax = param_estimates$nmax[b180_6_3_index], 
                                  mumax = (param_estimates$mumax[b180_6_3_index]*2.303),
                                  lag = param_estimates$lag[b180_6_3_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[73] <-"S12-0180"
fitted_param$pop[73] <-"apc"
fitted_param$model_type[73] <-"buchanan"
fitted_param$n0[73] <- summary(b180_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[73] <- summary(b180_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[73] <- summary(b180_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[73] <- summary(b180_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[73] <- "6"
fitted_param$aic[73] <- AIC(b180_buchanan_nls_lm_6_3)
fitted_param$bic[73] <- BIC(b180_buchanan_nls_lm_6_3)
fitted_param$fit[73] <- "yes"
fitted_param$method[73] <- "nls_lm"
fitted_param$rep[73] <- "3"

#buchanan no lag
b180_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b180_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b180_6_3_index], 
                                         log10nmax = param_estimates$nmax[b180_6_3_index], 
                                         mumax = (param_estimates$mumax[b180_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[74] <-"S12-0180"
fitted_param$pop[74] <-"apc"
fitted_param$model_type[74] <-"buchanan_no_lag"
fitted_param$n0[74] <- summary(b180_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[74] <- summary(b180_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[74] <- summary(b180_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[74] <- "6"
fitted_param$aic[74] <- AIC(b180_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[74] <- BIC(b180_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[74] <- "yes"
fitted_param$method[74] <- "nls_lm"
fitted_param$rep[74] <- "3"

#gompertz
b180_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b180_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b180_6_3_index], 
                                  log10nmax = param_estimates$nmax[b180_6_3_index], 
                                  mumax = (param_estimates$mumax[b180_6_3_index]*2.303),
                                  lag = param_estimates$lag[b180_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[75] <-"S12-0180"
fitted_param$pop[75] <-"apc"
fitted_param$model_type[75] <-"gompertz"
fitted_param$n0[75] <- summary(b180_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[75] <- summary(b180_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[75] <- summary(b180_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[75] <- summary(b180_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[75] <- "6"
fitted_param$aic[75] <- AIC(b180_gompertz_nls_lm_6_3)
fitted_param$bic[75] <- BIC(b180_gompertz_nls_lm_6_3)
fitted_param$fit[75] <- "yes"
fitted_param$method[75] <- "nls_lm"
fitted_param$rep[75] <- "3"

#End of S12-0180 model fitting 

##------------------------S12-0184 Model Fitting, 6C----------------------------

#Rep 1
b184_6_1_index <- which(param_estimates$isolate == "S12-0184" & param_estimates$temp == "6" & param_estimates$rep == "1")

#baranyi
b184_baranyi_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b184_6_1, 
                             start=list(
                               log10n0 = param_estimates$n0[b184_6_1_index], 
                               log10nmax = param_estimates$nmax[b184_6_1_index], 
                               mumax = (param_estimates$mumax[b184_6_1_index]*2.303), 
                               lag = param_estimates$lag[b184_6_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[76] <-"S12-0184"
fitted_param$pop[76] <-"apc"
fitted_param$model_type[76] <-"baranyi"
fitted_param$n0[76] <- summary(b184_baranyi_nls_lm_6_1)$coefficient[1]
fitted_param$lag[76] <- summary(b184_baranyi_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[76] <- summary(b184_baranyi_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[76] <- summary(b184_baranyi_nls_lm_6_1)$coefficient[2]
fitted_param$temp[76] <- "6"
fitted_param$aic[76] <- AIC(b184_baranyi_nls_lm_6_1)
fitted_param$bic[76] <- BIC(b184_baranyi_nls_lm_6_1)
fitted_param$fit[76] <- "yes"
fitted_param$method[76] <- "nls_lm"
fitted_param$rep[76] <- "1"

#baranyi no lag
b184_baranyi_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b184_6_1, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b184_6_1_index], 
                                      log10nmax = param_estimates$nmax[b184_6_1_index], 
                                      mumax = (param_estimates$mumax[b184_6_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[77] <-"S12-0184"
fitted_param$pop[77] <-"apc"
fitted_param$model_type[77] <-"baranyi_no_lag"
fitted_param$n0[77] <- summary(b184_baranyi_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[77] <- summary(b184_baranyi_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[77] <- summary(b184_baranyi_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[77] <- "6"
fitted_param$aic[77] <- AIC(b184_baranyi_no_lag_nls_lm_6_1)
fitted_param$bic[77] <- BIC(b184_baranyi_no_lag_nls_lm_6_1)
fitted_param$fit[77] <- "yes"
fitted_param$method[77] <- "nls_lm"
fitted_param$rep[77] <- "1"

#buchanan
b184_buchanan_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b184_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b184_6_1_index], 
                                log10nmax = param_estimates$nmax[b184_6_1_index], 
                                mumax = (param_estimates$mumax[b184_6_1_index]*2.303), 
                                lag = param_estimates$lag[b184_6_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[78] <-"S12-0184"
fitted_param$pop[78] <-"apc"
fitted_param$model_type[78] <-"buchanan"
fitted_param$n0[78] <- summary(b184_buchanan_nls_lm_6_1)$coefficient[1]
fitted_param$lag[78] <- summary(b184_buchanan_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[78] <- summary(b184_buchanan_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[78] <- summary(b184_buchanan_nls_lm_6_1)$coefficient[2]
fitted_param$temp[78] <- "6"
fitted_param$aic[78] <- AIC(b184_buchanan_nls_lm_6_1)
fitted_param$bic[78] <- BIC(b184_buchanan_nls_lm_6_1)
fitted_param$fit[78] <- "yes"
fitted_param$method[78] <- "nls_lm"
fitted_param$rep[78] <- "1"

#buchanan no lag
b184_buchanan_no_lag_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b184_6_1, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b184_6_1_index], 
                                       log10nmax = param_estimates$nmax[b184_6_1_index], 
                                       mumax = (param_estimates$mumax[b184_6_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[79] <-"S12-0184"
fitted_param$pop[79] <-"apc"
fitted_param$model_type[79] <-"buchanan_no_lag"
fitted_param$n0[79] <- summary(b184_buchanan_no_lag_nls_lm_6_1)$coefficient[1]
fitted_param$mumax[79] <- summary(b184_buchanan_no_lag_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[79] <- summary(b184_buchanan_no_lag_nls_lm_6_1)$coefficient[2]
fitted_param$temp[79] <- "6"
fitted_param$aic[79] <- AIC(b184_buchanan_no_lag_nls_lm_6_1)
fitted_param$bic[79] <- BIC(b184_buchanan_no_lag_nls_lm_6_1)
fitted_param$fit[79] <- "yes"
fitted_param$method[79] <- "nls_lm"
fitted_param$rep[79] <- "1"

#gompertz
b184_gompertz_nls_lm_6_1 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b184_6_1, 
                              start=list(
                                log10n0 = param_estimates$n0[b184_6_1_index], 
                                log10nmax = param_estimates$nmax[b184_6_1_index], 
                                mumax = (param_estimates$mumax[b184_6_1_index]*2.303),
                                lag = param_estimates$lag[b184_6_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[80] <-"S12-0184"
fitted_param$pop[80] <-"apc"
fitted_param$model_type[80] <-"gompertz"
fitted_param$n0[80] <- summary(b184_gompertz_nls_lm_6_1)$coefficient[1]
fitted_param$lag[80] <- summary(b184_gompertz_nls_lm_6_1)$coefficient[4]
fitted_param$mumax[80] <- summary(b184_gompertz_nls_lm_6_1)$coefficient[3]
fitted_param$nmax[80] <- summary(b184_gompertz_nls_lm_6_1)$coefficient[2]
fitted_param$temp[80] <- "6"
fitted_param$aic[80] <- AIC(b184_gompertz_nls_lm_6_1)
fitted_param$bic[80] <- BIC(b184_gompertz_nls_lm_6_1)
fitted_param$fit[80] <- "yes"
fitted_param$method[80] <- "nls_lm"
fitted_param$rep[80] <- "1"

#Rep 2
b184_6_2_index <- which(param_estimates$isolate == "S12-0184" & param_estimates$temp == "6" & param_estimates$rep == "2")

#baranyi
b184_baranyi_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b184_6_2, 
                               start=list(
                                 log10n0 = param_estimates$n0[b184_6_2_index], 
                                 log10nmax = param_estimates$nmax[b184_6_2_index], 
                                 mumax = (param_estimates$mumax[b184_6_2_index]*2.303), 
                                 lag = param_estimates$lag[b184_6_2_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[81] <-"S12-0184"
fitted_param$pop[81] <-"apc"
fitted_param$model_type[81] <-"baranyi"
fitted_param$n0[81] <- summary(b184_baranyi_nls_lm_6_2)$coefficient[1]
fitted_param$lag[81] <- summary(b184_baranyi_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[81] <- summary(b184_baranyi_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[81] <- summary(b184_baranyi_nls_lm_6_2)$coefficient[2]
fitted_param$temp[81] <- "6"
fitted_param$aic[81] <- AIC(b184_baranyi_nls_lm_6_2)
fitted_param$bic[81] <- BIC(b184_baranyi_nls_lm_6_2)
fitted_param$fit[81] <- "yes"
fitted_param$method[81] <- "nls_lm"
fitted_param$rep[81] <- "2"

#baranyi no lag
b184_baranyi_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b184_6_2, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b184_6_2_index], 
                                        log10nmax = param_estimates$nmax[b184_6_2_index], 
                                        mumax = (param_estimates$mumax[b184_6_2_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[82] <-"S12-0184"
fitted_param$pop[82] <-"apc"
fitted_param$model_type[82] <-"baranyi_no_lag"
fitted_param$n0[82] <- summary(b184_baranyi_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[82] <- summary(b184_baranyi_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[82] <- summary(b184_baranyi_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[82] <- "6"
fitted_param$aic[82] <- AIC(b184_baranyi_no_lag_nls_lm_6_2)
fitted_param$bic[82] <- BIC(b184_baranyi_no_lag_nls_lm_6_2)
fitted_param$fit[82] <- "yes"
fitted_param$method[82] <- "nls_lm"
fitted_param$rep[82] <- "2"

#buchanan
b184_buchanan_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b184_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b184_6_2_index], 
                                  log10nmax = param_estimates$nmax[b184_6_2_index], 
                                  mumax = (param_estimates$mumax[b184_6_2_index]*2.303), 
                                  lag = param_estimates$lag[b184_6_2_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[83] <-"S12-0184"
fitted_param$pop[83] <-"apc"
fitted_param$model_type[83] <-"buchanan"
fitted_param$n0[83] <- summary(b184_buchanan_nls_lm_6_2)$coefficient[1]
fitted_param$lag[83] <- summary(b184_buchanan_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[83] <- summary(b184_buchanan_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[83] <- summary(b184_buchanan_nls_lm_6_2)$coefficient[2]
fitted_param$temp[83] <- "6"
fitted_param$aic[83] <- AIC(b184_buchanan_nls_lm_6_2)
fitted_param$bic[83] <- BIC(b184_buchanan_nls_lm_6_2)
fitted_param$fit[83] <- "yes"
fitted_param$method[83] <- "nls_lm"
fitted_param$rep[83] <- "2"

#buchanan no lag
b184_buchanan_no_lag_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b184_6_2, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b184_6_2_index], 
                                         log10nmax = param_estimates$nmax[b184_6_2_index], 
                                         mumax = (param_estimates$mumax[b184_6_2_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[84] <-"S12-0184"
fitted_param$pop[84] <-"apc"
fitted_param$model_type[84] <-"buchanan_no_lag"
fitted_param$n0[84] <- summary(b184_buchanan_no_lag_nls_lm_6_2)$coefficient[1]
fitted_param$mumax[84] <- summary(b184_buchanan_no_lag_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[84] <- summary(b184_buchanan_no_lag_nls_lm_6_2)$coefficient[2]
fitted_param$temp[84] <- "6"
fitted_param$aic[84] <- AIC(b184_buchanan_no_lag_nls_lm_6_2)
fitted_param$bic[84] <- BIC(b184_buchanan_no_lag_nls_lm_6_2)
fitted_param$fit[84] <- "yes"
fitted_param$method[84] <- "nls_lm"
fitted_param$rep[84] <- "2"

#gompertz
b184_gompertz_nls_lm_6_2 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b184_6_2, 
                                start=list(
                                  log10n0 = param_estimates$n0[b184_6_2_index], 
                                  log10nmax = param_estimates$nmax[b184_6_2_index], 
                                  mumax = (param_estimates$mumax[b184_6_2_index]*2.303),
                                  lag = param_estimates$lag[b184_6_2_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[85] <-"S12-0184"
fitted_param$pop[85] <-"apc"
fitted_param$model_type[85] <-"gompertz"
fitted_param$n0[85] <- summary(b184_gompertz_nls_lm_6_2)$coefficient[1]
fitted_param$lag[85] <- summary(b184_gompertz_nls_lm_6_2)$coefficient[4]
fitted_param$mumax[85] <- summary(b184_gompertz_nls_lm_6_2)$coefficient[3]
fitted_param$nmax[85] <- summary(b184_gompertz_nls_lm_6_2)$coefficient[2]
fitted_param$temp[85] <- "6"
fitted_param$aic[85] <- AIC(b184_gompertz_nls_lm_6_2)
fitted_param$bic[85] <- BIC(b184_gompertz_nls_lm_6_2)
fitted_param$fit[85] <- "yes"
fitted_param$method[85] <- "nls_lm"
fitted_param$rep[85] <- "2"

#Rep 3
b184_6_3_index <- which(param_estimates$isolate == "S12-0184" & param_estimates$temp == "6" & param_estimates$rep == "3")

#baranyi
b184_baranyi_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                 baranyi(day, log10n0, log10nmax, mumax, lag), 
                               data = b184_6_3, 
                               start=list(
                                 log10n0 = param_estimates$n0[b184_6_3_index], 
                                 log10nmax = param_estimates$nmax[b184_6_3_index], 
                                 mumax = (param_estimates$mumax[b184_6_3_index]*2.303), 
                                 lag = param_estimates$lag[b184_6_3_index]),
                               lower = c(0, 0, 0, 0))

fitted_param$isolate[86] <-"S12-0184"
fitted_param$pop[86] <-"apc"
fitted_param$model_type[86] <-"baranyi"
fitted_param$n0[86] <- summary(b184_baranyi_nls_lm_6_3)$coefficient[1]
fitted_param$lag[86] <- summary(b184_baranyi_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[86] <- summary(b184_baranyi_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[86] <- summary(b184_baranyi_nls_lm_6_3)$coefficient[2]
fitted_param$temp[86] <- "6"
fitted_param$aic[86] <- AIC(b184_baranyi_nls_lm_6_3)
fitted_param$bic[86] <- BIC(b184_baranyi_nls_lm_6_3)
fitted_param$fit[86] <- "no"
fitted_param$method[86] <- "nls_lm"
fitted_param$rep[86] <- "3"

#baranyi no lag
b184_baranyi_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                        baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                      data = b184_6_3, 
                                      start=list(
                                        log10n0 = param_estimates$n0[b184_6_3_index], 
                                        log10nmax = param_estimates$nmax[b184_6_3_index], 
                                        mumax = (param_estimates$mumax[b184_6_3_index]*2.303)),
                                      lower = c(0, 0, 0))

fitted_param$isolate[87] <-"S12-0184"
fitted_param$pop[87] <-"apc"
fitted_param$model_type[87] <-"baranyi_no_lag"
fitted_param$n0[87] <- summary(b184_baranyi_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[87] <- summary(b184_baranyi_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[87] <- summary(b184_baranyi_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[87] <- "6"
fitted_param$aic[87] <- AIC(b184_baranyi_no_lag_nls_lm_6_3)
fitted_param$bic[87] <- BIC(b184_baranyi_no_lag_nls_lm_6_3)
fitted_param$fit[87] <- "no"
fitted_param$method[87] <- "nls_lm"
fitted_param$rep[87] <- "3"

#buchanan
b184_buchanan_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  buchanan(day, log10n0, log10nmax, mumax, lag), 
                                data = b184_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b184_6_3_index], 
                                  log10nmax = param_estimates$nmax[b184_6_3_index], 
                                  mumax = (param_estimates$mumax[b184_6_3_index]*2.303), 
                                  lag = param_estimates$lag[b184_6_3_index]), 
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[88] <-"S12-0184"
fitted_param$pop[88] <-"apc"
fitted_param$model_type[88] <-"buchanan"
fitted_param$n0[88] <- summary(b184_buchanan_nls_lm_6_3)$coefficient[1]
fitted_param$lag[88] <- summary(b184_buchanan_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[88] <- summary(b184_buchanan_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[88] <- summary(b184_buchanan_nls_lm_6_3)$coefficient[2]
fitted_param$temp[88] <- "6"
fitted_param$aic[88] <- AIC(b184_buchanan_nls_lm_6_3)
fitted_param$bic[88] <- BIC(b184_buchanan_nls_lm_6_3)
fitted_param$fit[88] <- "yes"
fitted_param$method[88] <- "nls_lm"
fitted_param$rep[88] <- "3"

#buchanan no lag
b184_buchanan_no_lag_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                         buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                       data = b184_6_3, 
                                       start=list(
                                         log10n0 = param_estimates$n0[b184_6_3_index], 
                                         log10nmax = param_estimates$nmax[b184_6_3_index], 
                                         mumax = (param_estimates$mumax[b184_6_3_index]*2.303)),
                                       lower = c(0, 0, 0))

fitted_param$isolate[89] <-"S12-0184"
fitted_param$pop[89] <-"apc"
fitted_param$model_type[89] <-"buchanan_no_lag"
fitted_param$n0[89] <- summary(b184_buchanan_no_lag_nls_lm_6_3)$coefficient[1]
fitted_param$mumax[89] <- summary(b184_buchanan_no_lag_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[89] <- summary(b184_buchanan_no_lag_nls_lm_6_3)$coefficient[2]
fitted_param$temp[89] <- "6"
fitted_param$aic[89] <- AIC(b184_buchanan_no_lag_nls_lm_6_3)
fitted_param$bic[89] <- BIC(b184_buchanan_no_lag_nls_lm_6_3)
fitted_param$fit[89] <- "yes"
fitted_param$method[89] <- "nls_lm"
fitted_param$rep[89] <- "3"

#gompertz
b184_gompertz_nls_lm_6_3 <- nlsLM(log_average_wrangled_conc ~ 
                                  gompertz(day, log10n0, log10nmax, mumax, lag), 
                                data = b184_6_3, 
                                start=list(
                                  log10n0 = param_estimates$n0[b184_6_3_index], 
                                  log10nmax = param_estimates$nmax[b184_6_3_index], 
                                  mumax = (param_estimates$mumax[b184_6_3_index]*2.303),
                                  lag = param_estimates$lag[b184_6_3_index]),
                                lower = c(0, 0, 0, 0))

fitted_param$isolate[90] <-"S12-0184"
fitted_param$pop[90] <-"apc"
fitted_param$model_type[90] <-"gompertz"
fitted_param$n0[90] <- summary(b184_gompertz_nls_lm_6_3)$coefficient[1]
fitted_param$lag[90] <- summary(b184_gompertz_nls_lm_6_3)$coefficient[4]
fitted_param$mumax[90] <- summary(b184_gompertz_nls_lm_6_3)$coefficient[3]
fitted_param$nmax[90] <- summary(b184_gompertz_nls_lm_6_3)$coefficient[2]
fitted_param$temp[90] <- "6"
fitted_param$aic[90] <- AIC(b184_gompertz_nls_lm_6_3)
fitted_param$bic[90] <- BIC(b184_gompertz_nls_lm_6_3)
fitted_param$fit[90] <- "yes"
fitted_param$method[90] <- "nls_lm"
fitted_param$rep[90] <- "3"

#End of S12-0184 model fitting 

##---------------------------S12-0116 Model Fitting, 10C-------------------------
b116_10_1_index <- which(param_estimates$isolate == "S12-0116" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b116_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b116_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b116_10_1_index], 
                               log10nmax = param_estimates$nmax[b116_10_1_index], 
                               mumax = (param_estimates$mumax[b116_10_1_index]*2.303), 
                               lag = param_estimates$lag[b116_10_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[91] <-"S12-0116"
fitted_param$pop[91] <-"apc"
fitted_param$model_type[91] <-"baranyi"
fitted_param$n0[91] <- summary(b116_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[91] <- summary(b116_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[91] <- summary(b116_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[91] <- summary(b116_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[91] <- "10"
fitted_param$aic[91] <- AIC(b116_baranyi_nls_lm_10)
fitted_param$bic[91] <- BIC(b116_baranyi_nls_lm_10)
fitted_param$fit[91] <- "yes"
fitted_param$method[91] <- "nls_lm"
fitted_param$rep[91] <- "1"

#baranyi no lag
b116_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b116_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b116_10_1_index], 
                                      log10nmax = param_estimates$nmax[b116_10_1_index], 
                                      mumax = (param_estimates$mumax[b116_10_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[92] <-"S12-0116"
fitted_param$pop[92] <-"apc"
fitted_param$model_type[92] <-"baranyi_no_lag"
fitted_param$n0[92] <- summary(b116_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[92] <- summary(b116_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[92] <- summary(b116_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[92] <- "10"
fitted_param$aic[92] <- AIC(b116_baranyi_no_lag_nls_lm_10)
fitted_param$bic[92] <- BIC(b116_baranyi_no_lag_nls_lm_10)
fitted_param$fit[92] <- "yes"
fitted_param$method[92] <- "nls_lm"
fitted_param$rep[92] <- "1"

#buchanan
b116_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b116_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b116_10_1_index], 
                                log10nmax = param_estimates$nmax[b116_10_1_index], 
                                mumax = (param_estimates$mumax[b116_10_1_index]*2.303), 
                                lag = param_estimates$lag[b116_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[93] <-"S12-0116"
fitted_param$pop[93] <-"apc"
fitted_param$model_type[93] <-"buchanan"
fitted_param$n0[93] <- summary(b116_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[93] <- summary(b116_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[93] <- summary(b116_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[93] <- summary(b116_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[93] <- "10"
fitted_param$aic[93] <- AIC(b116_buchanan_nls_lm_10)
fitted_param$bic[93] <- BIC(b116_buchanan_nls_lm_10)
fitted_param$fit[93] <- "yes"
fitted_param$method[93] <- "nls_lm"
fitted_param$rep[93] <- "1"

#buchanan no lag
b116_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b116_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b116_10_1_index], 
                                       log10nmax = param_estimates$nmax[b116_10_1_index], 
                                       mumax = (param_estimates$mumax[b116_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[94] <-"S12-0116"
fitted_param$pop[94] <-"apc"
fitted_param$model_type[94] <-"buchanan_no_lag"
fitted_param$n0[94] <- summary(b116_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[94] <- summary(b116_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[94] <- summary(b116_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[94] <- "10"
fitted_param$aic[94] <- AIC(b116_buchanan_no_lag_nls_lm_10)
fitted_param$bic[94] <- BIC(b116_buchanan_no_lag_nls_lm_10)
fitted_param$fit[94] <- "yes"
fitted_param$method[94] <- "nls_lm"
fitted_param$rep[94] <- "1"

#gompertz
b116_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b116_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b116_10_1_index], 
                                log10nmax = param_estimates$nmax[b116_10_1_index], 
                                mumax = (param_estimates$mumax[b116_10_1_index]*2.303),
                                lag = param_estimates$lag[b116_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[95] <-"S12-0116"
fitted_param$pop[95] <-"apc"
fitted_param$model_type[95] <-"gompertz"
fitted_param$n0[95] <- summary(b116_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[95] <- summary(b116_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[95] <- summary(b116_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[95] <- summary(b116_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[95] <- "10"
fitted_param$aic[95] <- AIC(b116_gompertz_nls_lm_10)
fitted_param$bic[95] <- BIC(b116_gompertz_nls_lm_10)
fitted_param$fit[95] <- "yes"
fitted_param$method[95] <- "nls_lm"
fitted_param$rep[95] <- "1"

#End of S12-0116 model fitting 

##-------------------------S12-0132 Model Fitting, 10C---------------------------
b132_10_1_index <- which(param_estimates$isolate == "S12-0132" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b132_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b132_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b132_10_1_index], 
                               log10nmax = param_estimates$nmax[b132_10_1_index], 
                               mumax = (param_estimates$mumax[b132_10_1_index]*2.303), 
                               lag = param_estimates$lag[b132_10_1_index]), 
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[96] <-"S12-0132"
fitted_param$pop[96] <-"apc"
fitted_param$model_type[96] <-"baranyi"
fitted_param$n0[96] <- summary(b132_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[96] <- summary(b132_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[96] <- summary(b132_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[96] <- summary(b132_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[96] <- "10"
fitted_param$aic[96] <- AIC(b132_baranyi_nls_lm_10)
fitted_param$bic[96] <- BIC(b132_baranyi_nls_lm_10)
fitted_param$fit[96] <- "yes"
fitted_param$method[96] <- "nls_lm"
fitted_param$rep[96] <- "1"

#baranyi no lag
b132_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b132_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b132_10_1_index], 
                                      log10nmax = param_estimates$nmax[b132_10_1_index], 
                                      mumax = (param_estimates$mumax[b132_10_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[97] <-"S12-0132"
fitted_param$pop[97] <-"apc"
fitted_param$model_type[97] <-"baranyi_no_lag"
fitted_param$n0[97] <- summary(b132_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[97] <- summary(b132_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[97] <- summary(b132_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[97] <- "10"
fitted_param$aic[97] <- AIC(b132_baranyi_no_lag_nls_lm_10)
fitted_param$bic[97] <- BIC(b132_baranyi_no_lag_nls_lm_10)
fitted_param$fit[97] <- "yes"
fitted_param$method[97] <- "nls_lm"
fitted_param$rep[97] <- "1"

#buchanan
b132_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b132_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b132_10_1_index], 
                                log10nmax = param_estimates$nmax[b132_10_1_index], 
                                mumax = (param_estimates$mumax[b132_10_1_index]*2.303), 
                                lag = param_estimates$lag[b132_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[98] <-"S12-0132"
fitted_param$pop[98] <-"apc"
fitted_param$model_type[98] <-"buchanan"
fitted_param$n0[98] <- summary(b132_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[98] <- summary(b132_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[98] <- summary(b132_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[98] <- summary(b132_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[98] <- "10"
fitted_param$aic[98] <- AIC(b132_buchanan_nls_lm_10)
fitted_param$bic[98] <- BIC(b132_buchanan_nls_lm_10)
fitted_param$fit[98] <- "yes"
fitted_param$method[98] <- "nls_lm"
fitted_param$rep[98] <- "1"

#buchanan no lag
b132_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b132_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b132_10_1_index], 
                                       log10nmax = param_estimates$nmax[b132_10_1_index], 
                                       mumax = (param_estimates$mumax[b132_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[99] <-"S12-0132"
fitted_param$pop[99] <-"apc"
fitted_param$model_type[99] <-"buchanan_no_lag"
fitted_param$n0[99] <- summary(b132_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[99] <- summary(b132_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[99] <- summary(b132_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[99] <- "10"
fitted_param$aic[99] <- AIC(b132_buchanan_no_lag_nls_lm_10)
fitted_param$bic[99] <- BIC(b132_buchanan_no_lag_nls_lm_10)
fitted_param$fit[99] <- "yes"
fitted_param$method[99] <- "nls_lm"
fitted_param$rep[99] <- "1"

#gompertz
b132_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b132_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b132_10_1_index], 
                                log10nmax = param_estimates$nmax[b132_10_1_index], 
                                mumax = (param_estimates$mumax[b132_10_1_index]*2.303),
                                lag = param_estimates$lag[b132_10_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[100] <-"S12-0132"
fitted_param$pop[100] <-"apc"
fitted_param$model_type[100] <-"gompertz"
fitted_param$n0[100] <- summary(b132_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[100] <- summary(b132_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[100] <- summary(b132_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[100] <- summary(b132_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[100] <- "10"
fitted_param$aic[100] <- AIC(b132_gompertz_nls_lm_10)
fitted_param$bic[100] <- BIC(b132_gompertz_nls_lm_10)
fitted_param$fit[100] <- "yes"
fitted_param$method[100] <- "nls_lm"
fitted_param$rep[100] <- "1"

#End of S12-0132 model fitting 

##------------------------S12-0141 Model Fitting, 10C----------------------------
b141_10_1_index <- which(param_estimates$isolate == "S12-0141" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b141_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b141_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b141_10_1_index], 
                               log10nmax = param_estimates$nmax[b141_10_1_index], 
                               mumax = (param_estimates$mumax[b141_10_1_index]*2.303), 
                               lag = param_estimates$lag[b141_10_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[101] <-"S12-0141"
fitted_param$pop[101] <-"apc"
fitted_param$model_type[101] <-"baranyi"
fitted_param$n0[101] <- summary(b141_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[101] <- summary(b141_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[101] <- summary(b141_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[101] <- summary(b141_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[101] <- "10"
fitted_param$aic[101] <- AIC(b141_baranyi_nls_lm_10)
fitted_param$bic[101] <- BIC(b141_baranyi_nls_lm_10)
fitted_param$fit[101] <- "yes"
fitted_param$method[101] <- "nls_lm"
fitted_param$rep[101] <- "1"

#baranyi no lag
b141_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b141_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b141_10_1_index], 
                                      log10nmax = param_estimates$nmax[b141_10_1_index], 
                                      mumax = (param_estimates$mumax[b141_10_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[102] <-"S12-0141"
fitted_param$pop[102] <-"apc"
fitted_param$model_type[102] <-"baranyi_no_lag"
fitted_param$n0[102] <- summary(b141_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[102] <- summary(b141_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[102] <- summary(b141_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[102] <- "10"
fitted_param$aic[102] <- AIC(b141_baranyi_no_lag_nls_lm_10)
fitted_param$bic[102] <- BIC(b141_baranyi_no_lag_nls_lm_10)
fitted_param$fit[102] <- "yes"
fitted_param$method[102] <- "nls_lm"
fitted_param$rep[102] <- "1"

#buchanan
b141_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b141_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b141_10_1_index], 
                                log10nmax = param_estimates$nmax[b141_10_1_index], 
                                mumax = (param_estimates$mumax[b141_10_1_index]*2.303), 
                                lag = param_estimates$lag[b141_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[103] <-"S12-0141"
fitted_param$pop[103] <-"apc"
fitted_param$model_type[103] <-"buchanan"
fitted_param$n0[103] <- summary(b141_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[103] <- summary(b141_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[103] <- summary(b141_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[103] <- summary(b141_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[103] <- "10"
fitted_param$aic[103] <- AIC(b141_buchanan_nls_lm_10)
fitted_param$bic[103] <- BIC(b141_buchanan_nls_lm_10)
fitted_param$fit[103] <- "yes"
fitted_param$method[103] <- "nls_lm"
fitted_param$rep[103] <- "1"

#buchanan no lag
b141_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b141_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b141_10_1_index], 
                                       log10nmax = param_estimates$nmax[b141_10_1_index], 
                                       mumax = (param_estimates$mumax[b141_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[104] <-"S12-0141"
fitted_param$pop[104] <-"apc"
fitted_param$model_type[104] <-"buchanan_no_lag"
fitted_param$n0[104] <- summary(b141_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[104] <- summary(b141_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[104] <- summary(b141_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[104] <- "10"
fitted_param$aic[104] <- AIC(b141_buchanan_no_lag_nls_lm_10)
fitted_param$bic[104] <- BIC(b141_buchanan_no_lag_nls_lm_10)
fitted_param$fit[104] <- "yes"
fitted_param$method[104] <- "nls_lm"
fitted_param$rep[104] <- "1"

#gompertz
b141_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b141_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b141_10_1_index], 
                                log10nmax = param_estimates$nmax[b141_10_1_index], 
                                mumax = (param_estimates$mumax[b141_10_1_index]*2.303),
                                lag = param_estimates$lag[b141_10_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[105] <-"S12-0141"
fitted_param$pop[105] <-"apc"
fitted_param$model_type[105] <-"gompertz"
fitted_param$n0[105] <- summary(b141_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[105] <- summary(b141_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[105] <- summary(b141_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[105] <- summary(b141_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[105] <- "10"
fitted_param$aic[105] <- AIC(b141_gompertz_nls_lm_10)
fitted_param$bic[105] <- BIC(b141_gompertz_nls_lm_10)
fitted_param$fit[105] <- "yes"
fitted_param$method[105] <- "nls_lm"
fitted_param$rep[105] <- "1"

#End of S12-0141 model fitting 

##------------------------S12-0166 Model Fitting, 10C---------------------------
b166_10_1_index <- which(param_estimates$isolate == "S12-0166" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b166_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b166_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b166_10_1_index], 
                               log10nmax = param_estimates$nmax[b166_10_1_index], 
                               mumax = (param_estimates$mumax[b166_10_1_index]*2.303), 
                               lag = param_estimates$lag[b166_10_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[106] <-"S12-0166"
fitted_param$pop[106] <-"apc"
fitted_param$model_type[106] <-"baranyi"
fitted_param$n0[106] <- summary(b166_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[106] <- summary(b166_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[106] <- summary(b166_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[106] <- summary(b166_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[106] <- "10"
fitted_param$aic[106] <- AIC(b166_baranyi_nls_lm_10)
fitted_param$bic[106] <- BIC(b166_baranyi_nls_lm_10)
fitted_param$fit[106] <- "yes"
fitted_param$method[106] <- "nls_lm"
fitted_param$rep[106] <- "1"

#baranyi no lag
b166_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b166_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b166_10_1_index], 
                                      log10nmax = param_estimates$nmax[b166_10_1_index], 
                                      mumax = (param_estimates$mumax[b166_10_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[107] <-"S12-0166"
fitted_param$pop[107] <-"apc"
fitted_param$model_type[107] <-"baranyi_no_lag"
fitted_param$n0[107] <- summary(b166_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[107] <- summary(b166_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[107] <- summary(b166_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[107] <- "10"
fitted_param$aic[107] <- AIC(b166_baranyi_no_lag_nls_lm_10)
fitted_param$bic[107] <- BIC(b166_baranyi_no_lag_nls_lm_10)
fitted_param$fit[107] <- "yes"
fitted_param$method[107] <- "nls_lm"
fitted_param$rep[107] <- "1"

#buchanan
b166_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b166_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b166_10_1_index], 
                                log10nmax = param_estimates$nmax[b166_10_1_index], 
                                mumax = (param_estimates$mumax[b166_10_1_index]*2.303), 
                                lag = param_estimates$lag[b166_10_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[108] <-"S12-0166"
fitted_param$pop[108] <-"apc"
fitted_param$model_type[108] <-"buchanan"
fitted_param$n0[108] <- summary(b166_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[108] <- summary(b166_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[108] <- summary(b166_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[108] <- summary(b166_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[108] <- "10"
fitted_param$aic[108] <- AIC(b166_buchanan_nls_lm_10)
fitted_param$bic[108] <- BIC(b166_buchanan_nls_lm_10)
fitted_param$fit[108] <- "yes"
fitted_param$method[108] <- "nls_lm"
fitted_param$rep[108] <- "1"

#buchanan no lag
b166_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b166_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b166_10_1_index], 
                                       log10nmax = param_estimates$nmax[b166_10_1_index], 
                                       mumax = (param_estimates$mumax[b166_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[109] <-"S12-0166"
fitted_param$pop[109] <-"apc"
fitted_param$model_type[109] <-"buchanan_no_lag"
fitted_param$n0[109] <- summary(b166_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[109] <- summary(b166_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[109] <- summary(b166_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[109] <- "10"
fitted_param$aic[109] <- AIC(b166_buchanan_no_lag_nls_lm_10)
fitted_param$bic[109] <- BIC(b166_buchanan_no_lag_nls_lm_10)
fitted_param$fit[109] <- "yes"
fitted_param$method[109] <- "nls_lm"
fitted_param$rep[109] <- "1"

#gompertz
b166_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b166_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b166_10_1_index], 
                                log10nmax = param_estimates$nmax[b166_10_1_index], 
                                mumax = (param_estimates$mumax[b166_10_1_index]*2.303),
                                lag = param_estimates$lag[b166_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[110] <-"S12-0166"
fitted_param$pop[110] <-"apc"
fitted_param$model_type[110] <-"gompertz"
fitted_param$n0[110] <- summary(b166_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[110] <- summary(b166_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[110] <- summary(b166_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[110] <- summary(b166_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[110] <- "10"
fitted_param$aic[110] <- AIC(b166_gompertz_nls_lm_10)
fitted_param$bic[110] <- BIC(b166_gompertz_nls_lm_10)
fitted_param$fit[110] <- "yes"
fitted_param$method[110] <- "nls_lm"
fitted_param$rep[110] <- "1"

#End of S12-0166 model fitting 

##------------------------S12-0180 Model Fitting, 10C---------------------------
b180_10_1_index <- which(param_estimates$isolate == "S12-0180" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b180_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b180_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b180_10_1_index], 
                               log10nmax = param_estimates$nmax[b180_10_1_index], 
                               mumax = (param_estimates$mumax[b180_10_1_index]*2.303), 
                               lag = param_estimates$lag[b180_10_1_index]),
                             lower = c(0, 0, 0, 0))

fitted_param$isolate[111] <-"S12-0180"
fitted_param$pop[111] <-"apc"
fitted_param$model_type[111] <-"baranyi"
fitted_param$n0[111] <- summary(b180_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[111] <- summary(b180_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[111] <- summary(b180_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[111] <- summary(b180_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[111] <- "10"
fitted_param$aic[111] <- AIC(b180_baranyi_nls_lm_10)
fitted_param$bic[111] <- BIC(b180_baranyi_nls_lm_10)
fitted_param$fit[111] <- "yes_maxiter_error"
fitted_param$method[111] <- "nls_lm"
fitted_param$rep[111] <- "1"

#baranyi no lag
b180_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b180_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b180_10_1_index], 
                                      log10nmax = param_estimates$nmax[b180_10_1_index], 
                                      mumax = param_estimates$mumax[b180_10_1_index]),
                                    lower = c(0, 0, 0))

fitted_param$isolate[112] <-"S12-0180"
fitted_param$pop[112] <-"apc"
fitted_param$model_type[112] <-"baranyi_no_lag"
fitted_param$n0[112] <- summary(b180_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[112] <- summary(b180_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[112] <- summary(b180_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[112] <- "10"
fitted_param$aic[112] <- AIC(b180_baranyi_no_lag_nls_lm_10)
fitted_param$bic[112] <- BIC(b180_baranyi_no_lag_nls_lm_10)
fitted_param$fit[112] <- "no"
fitted_param$method[112] <- "nls_lm"
fitted_param$rep[112] <- "1"

#buchanan
b180_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b180_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b180_10_1_index], 
                                log10nmax = param_estimates$nmax[b180_10_1_index], 
                                mumax = (param_estimates$mumax[b180_10_1_index]*2.303), 
                                lag = param_estimates$lag[b180_10_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[113] <-"S12-0180"
fitted_param$pop[113] <-"apc"
fitted_param$model_type[113] <-"buchanan"
fitted_param$n0[113] <- summary(b180_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[113] <- summary(b180_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[113] <- summary(b180_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[113] <- summary(b180_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[113] <- "10"
fitted_param$aic[113] <- AIC(b180_buchanan_nls_lm_10)
fitted_param$bic[113] <- BIC(b180_buchanan_nls_lm_10)
fitted_param$fit[113] <- "yes"
fitted_param$method[113] <- "nls_lm"
fitted_param$rep[113] <- "1"

#buchanan no lag
b180_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b180_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b180_10_1_index], 
                                       log10nmax = param_estimates$nmax[b180_10_1_index], 
                                       mumax = (param_estimates$mumax[b180_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[114] <-"S12-0180"
fitted_param$pop[114] <-"apc"
fitted_param$model_type[114] <-"buchanan_no_lag"
fitted_param$n0[114] <- summary(b180_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[114] <- summary(b180_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[114] <- summary(b180_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[114] <- "10"
fitted_param$aic[114] <- AIC(b180_buchanan_no_lag_nls_lm_10)
fitted_param$bic[114] <- BIC(b180_buchanan_no_lag_nls_lm_10)
fitted_param$fit[114] <- "yes"
fitted_param$method[114] <- "nls_lm"
fitted_param$rep[114] <- "1"

#gompertz
b180_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b180_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b180_10_1_index], 
                                log10nmax = param_estimates$nmax[b180_10_1_index], 
                                mumax = (param_estimates$mumax[b180_10_1_index]*2.303),
                                lag = param_estimates$lag[b180_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[115] <-"S12-0180"
fitted_param$pop[115] <-"apc"
fitted_param$model_type[115] <-"gompertz"
fitted_param$n0[115] <- summary(b180_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[115] <- summary(b180_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[115] <- summary(b180_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[115] <- summary(b180_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[115] <- "10"
fitted_param$aic[115] <- AIC(b180_gompertz_nls_lm_10)
fitted_param$bic[115] <- BIC(b180_gompertz_nls_lm_10)
fitted_param$fit[115] <- "yes"
fitted_param$method[115] <- "nls_lm"
fitted_param$rep[115] <- "1"

#End of S12-0180 model fitting 

##------------------------S12-0184 Model Fitting, 10C---------------------------
b184_10_1_index <- which(param_estimates$isolate == "S12-0184" & param_estimates$temp == "10" & param_estimates$rep == "1")

#baranyi
b184_baranyi_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                               baranyi(day, log10n0, log10nmax, mumax, lag), 
                             data = b184_10, 
                             start=list(
                               log10n0 = param_estimates$n0[b184_10_1_index], 
                               log10nmax = param_estimates$nmax[b184_10_1_index], 
                               mumax = (param_estimates$mumax[b184_10_1_index]*2.303), 
                               lag = param_estimates$lag[b184_10_1_index]),
                             lower = c(0, 0, 0, 0),
                             control = nls.lm.control(maxiter = 150))

fitted_param$isolate[116] <-"S12-0184"
fitted_param$pop[116] <-"apc"
fitted_param$model_type[116] <-"baranyi"
fitted_param$n0[116] <- summary(b184_baranyi_nls_lm_10)$coefficient[1]
fitted_param$lag[116] <- summary(b184_baranyi_nls_lm_10)$coefficient[4]
fitted_param$mumax[116] <- summary(b184_baranyi_nls_lm_10)$coefficient[3]
fitted_param$nmax[116] <- summary(b184_baranyi_nls_lm_10)$coefficient[2]
fitted_param$temp[116] <- "10"
fitted_param$aic[116] <- AIC(b184_baranyi_nls_lm_10)
fitted_param$bic[116] <- BIC(b184_baranyi_nls_lm_10)
fitted_param$fit[116] <- "no"
fitted_param$method[116] <- "nls_lm"
fitted_param$rep[116] <- "1"

#baranyi no lag
b184_baranyi_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                      baranyi_no_lag(day, log10n0, log10nmax, mumax), 
                                    data = b184_10, 
                                    start=list(
                                      log10n0 = param_estimates$n0[b184_10_1_index], 
                                      log10nmax = param_estimates$nmax[b184_10_1_index], 
                                      mumax = (param_estimates$mumax[b184_10_1_index]*2.303)),
                                    lower = c(0, 0, 0))

fitted_param$isolate[117] <-"S12-0184"
fitted_param$pop[117] <-"apc"
fitted_param$model_type[117] <-"baranyi_no_lag"
fitted_param$n0[117] <- summary(b184_baranyi_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[117] <- summary(b184_baranyi_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[117] <- summary(b184_baranyi_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[117] <- "10"
fitted_param$aic[117] <- AIC(b184_baranyi_no_lag_nls_lm_10)
fitted_param$bic[117] <- BIC(b184_baranyi_no_lag_nls_lm_10)
fitted_param$fit[117] <- "yes"
fitted_param$method[117] <- "nls_lm"
fitted_param$rep[117] <- "1"

#buchanan
b184_buchanan_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                buchanan(day, log10n0, log10nmax, mumax, lag), 
                              data = b184_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b184_10_1_index], 
                                log10nmax = param_estimates$nmax[b184_10_1_index], 
                                mumax = (param_estimates$mumax[b184_10_1_index]*2.303), 
                                lag = param_estimates$lag[b184_10_1_index]), 
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[118] <-"S12-0184"
fitted_param$pop[118] <-"apc"
fitted_param$model_type[118] <-"buchanan"
fitted_param$n0[118] <- summary(b184_buchanan_nls_lm_10)$coefficient[1]
fitted_param$lag[118] <- summary(b184_buchanan_nls_lm_10)$coefficient[4]
fitted_param$mumax[118] <- summary(b184_buchanan_nls_lm_10)$coefficient[3]
fitted_param$nmax[118] <- summary(b184_buchanan_nls_lm_10)$coefficient[2]
fitted_param$temp[118] <- "10"
fitted_param$aic[118] <- AIC(b184_buchanan_nls_lm_10)
fitted_param$bic[118] <- BIC(b184_buchanan_nls_lm_10)
fitted_param$fit[118] <- "yes"
fitted_param$method[118] <- "nls_lm"
fitted_param$rep[118] <- "1"

#buchanan no lag
b184_buchanan_no_lag_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                       buchanan_no_lag(day, log10n0, log10nmax, mumax), 
                                     data = b184_10, 
                                     start=list(
                                       log10n0 = param_estimates$n0[b184_10_1_index], 
                                       log10nmax = param_estimates$nmax[b184_10_1_index], 
                                       mumax = (param_estimates$mumax[b184_10_1_index]*2.303)),
                                     lower = c(0, 0, 0))

fitted_param$isolate[119] <-"S12-0184"
fitted_param$pop[119] <-"apc"
fitted_param$model_type[119] <-"buchanan_no_lag"
fitted_param$n0[119] <- summary(b184_buchanan_no_lag_nls_lm_10)$coefficient[1]
fitted_param$mumax[119] <- summary(b184_buchanan_no_lag_nls_lm_10)$coefficient[3]
fitted_param$nmax[119] <- summary(b184_buchanan_no_lag_nls_lm_10)$coefficient[2]
fitted_param$temp[119] <- "10"
fitted_param$aic[119] <- AIC(b184_buchanan_no_lag_nls_lm_10)
fitted_param$bic[119] <- BIC(b184_buchanan_no_lag_nls_lm_10)
fitted_param$fit[119] <- "yes"
fitted_param$method[119] <- "nls_lm"
fitted_param$rep[119] <- "1"

#gompertz
b184_gompertz_nls_lm_10 <- nlsLM(log_average_wrangled_conc ~ 
                                gompertz(day, log10n0, log10nmax, mumax, lag), 
                              data = b184_10, 
                              start=list(
                                log10n0 = param_estimates$n0[b184_10_1_index], 
                                log10nmax = param_estimates$nmax[b184_10_1_index], 
                                mumax = (param_estimates$mumax[b184_10_1_index]*2.303),
                                lag = param_estimates$lag[b184_10_1_index]),
                              lower = c(0, 0, 0, 0))

fitted_param$isolate[120] <-"S12-0184"
fitted_param$pop[120] <-"apc"
fitted_param$model_type[120] <-"gompertz"
fitted_param$n0[120] <- summary(b184_gompertz_nls_lm_10)$coefficient[1]
fitted_param$lag[120] <- summary(b184_gompertz_nls_lm_10)$coefficient[4]
fitted_param$mumax[120] <- summary(b184_gompertz_nls_lm_10)$coefficient[3]
fitted_param$nmax[120] <- summary(b184_gompertz_nls_lm_10)$coefficient[2]
fitted_param$temp[120] <- "10"
fitted_param$aic[120] <- AIC(b184_gompertz_nls_lm_10)
fitted_param$bic[120] <- BIC(b184_gompertz_nls_lm_10)
fitted_param$fit[120] <- "yes"
fitted_param$method[120] <- "nls_lm"
fitted_param$rep[120] <- "1"

#End of S12-0184 model fitting 
## ----------------------Selecting the model------------------------------------

buchanan_no_lag_aic <- fitted_param %>%
  filter(model_type == "buchanan_no_lag") %>%
  summarize(mean_aic = mean(aic),
            min_aic = min(aic),
            max_aic = max(aic))

gompertz_no_lag_aic <- fitted_param %>%
  filter(model_type == "gompertz") %>%
  summarize(mean_aic = mean(aic),
            min_aic = min(aic),
            max_aic = max(aic))

#Selecting the Buchanan-no-lag model, since all data could be fit to this model
#and because we are assuming there is no lag in this data
#All data could also be fit to the Gompertz model, but the fit (AIC) was higher, suggesting
#that is has a poorer fit than the Buchanan-no-lag model

final_model <- fitted_param %>%
  filter(model_type == "buchanan_no_lag") 

aic_by_model <- fitted_param %>%
  filter(model_type %in% c("buchanan_no_lag", "gompertz")) %>%
  group_by(model_type) %>%
  summarize(mean_aic = mean(aic)) 

set_1 <- c("S12-0116", "S12-0132", "S12-0141")
set_2 <- c("S12-0166", "S12-0180", "S12-0184")

param <- final_model %>%
  mutate(lot = case_when(
    (isolate %in% set_1 & rep == 1 & temp == 6) ~ 1,
    (isolate %in% set_1 & rep == 2 & temp == 6) ~ 2,
    (isolate %in% set_1 & rep == 3 & temp == 6) ~ 3,
    (isolate %in% set_2 & rep == 1 & temp == 6) ~ 4,
    (isolate %in% set_2 & rep == 2 & temp == 6) ~ 5,
    (isolate %in% set_2 & rep == 3 & temp == 6) ~ 6,
    (isolate %in% set_1 & temp == 10) ~ 7,
    (isolate %in% set_2 & temp == 10) ~ 8
  ))

param_ave <- param %>%
  group_by(lot) %>%
  mutate(mean_n0 = mean(n0), mean_lag = mean(lag), mean_mumax = mean(mumax), mean_nmax = mean(nmax)) %>%
  distinct(temp, pop, model_type, mean_n0, mean_lag, mean_mumax, mean_nmax, lot)

colnames(param_ave)[colnames(param_ave) %in% c("mean_n0", "mean_lag","mean_mumax", "mean_nmax")] <- 
  c("n0", "lag", "mumax", "nmax")

## -------------------------------Push data back onto server--------------------
# Save processed data to R Project folder 

#Push the processed data back to the R project 

#write.csv(final_model, "outputs/parameters/primary_model_parameters_apc_by_rep.csv", row.names = FALSE)
#write.csv(param_ave, "outputs/parameters/primary_model_parameters_apc_averaged.csv", row.names = FALSE)

# End of saving data into R Project folder 


