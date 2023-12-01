##  ----------------------------Title-------------------------------------------
#   Fitting the Ratkowsky square root model to the strain growth data (averaged growth parameters) 

##  --------------------------Description---------------------------------------
#   Project: CIDA Spinach 

#   Script description: Fitting a secondary model, to assess change in mu, based on primary growth parameters

##  --------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); 

##  --------------------------Data----------------------------------------------
#Read in data 
primary_param <- read.csv("outputs/parameters/primary_model_parameters_strain_averaged.csv", header = TRUE)

## ------------------Defining Secondary Model-----------------------------------

#Ratkowsky
ratkowsky <- function(temp, b, temp_min){
  sqrt_new_mu <- b*(temp - temp_min)
  new_mu <- sqrt_new_mu^2
  return(new_mu)
}

#End of defining secondary model

## ------------------Data Wrangling---------------------------------------------
#primary_param$temp <- weathermetrics::celsius.to.kelvin(primary_param$temp, 2)

b116 <- primary_param %>%
  filter(isolate == "S12-0116") %>%
  mutate(sqrt_mu = sqrt(mumax))

b132 <- primary_param %>%
  filter(isolate == "S12-0132") %>%
  mutate(sqrt_mu = sqrt(mumax)) 

b141 <- primary_param %>%
  filter(isolate == "S12-0141") %>%
  mutate(sqrt_mu = sqrt(mumax)) 

b166 <- primary_param %>%
  filter(isolate == "S12-0166") %>%
  mutate(sqrt_mu = sqrt(mumax)) 

b180 <- primary_param %>%
  filter(isolate == "S12-0180") %>%
  mutate(sqrt_mu = sqrt(mumax)) 

b184 <- primary_param %>%
  filter(isolate == "S12-0184") %>%
  mutate(sqrt_mu = sqrt(mumax)) 

#Making a dataframe to store the fitted model parameter estimates for the isolates 
no_of_isolates <- 6
secondary_param <- matrix(nrow = (no_of_isolates), ncol = 3)
secondary_param <- data.frame(secondary_param)
colnames(secondary_param) <- c("isolate", "b_mu", "temp_min_mu")

#End of data wrangling

## ------------------Model Fitting----------------------------------------------

#S12-0116
secondary_param[1, 1] <- "S12-0116"

b116_mod_mu <- lm(sqrt_mu ~ temp, data = b116)
secondary_param[1, 2] <- b116_mod_mu$coefficients[2]
secondary_param[1, 3] <- (-b116_mod_mu$coefficients[1])/b116_mod_mu$coefficients[2]


#S12-0132
secondary_param[2, 1] <- "S12-0132"

b132_mod_mu <- lm(sqrt_mu ~ temp, data = b132)
secondary_param[2, 2] <- b132_mod_mu$coefficients[2]
secondary_param[2, 3] <- (-b132_mod_mu$coefficients[1])/b132_mod_mu$coefficients[2]

#S12-0141
secondary_param[3, 1] <- "S12-0141"

b141_mod_mu <- lm(sqrt_mu ~ temp, data = b141)
secondary_param[3, 2] <- b141_mod_mu$coefficients[2]
secondary_param[3, 3] <- (-b141_mod_mu$coefficients[1])/b141_mod_mu$coefficients[2]

#S12-0166
secondary_param[4, 1] <- "S12-0166"

b166_mod_mu <- lm(sqrt_mu ~ temp, data = b166)
secondary_param[4, 2] <- b166_mod_mu$coefficients[2]
secondary_param[4, 3] <- (-b166_mod_mu$coefficients[1])/b166_mod_mu$coefficients[2]

#S12-0180
secondary_param[5, 1] <- "S12-0180"

b180_mod_mu <- lm(sqrt_mu ~ temp, data = b180)
secondary_param[5, 2] <- b180_mod_mu$coefficients[2]
secondary_param[5, 3] <- (-b180_mod_mu$coefficients[1])/b180_mod_mu$coefficients[2]

#S12-0184
secondary_param[6, 1] <- "S12-0184"

b184_mod_mu <- lm(sqrt_mu ~ temp, data = b184)
secondary_param[6, 2] <- b184_mod_mu$coefficients[2]
secondary_param[6, 3] <- (-b184_mod_mu$coefficients[1])/b184_mod_mu$coefficients[2]

## -------------------------------Push data back onto server--------------------
# Save data to R Project folder 

#Push the data back to the R project 

#Parameters
#write.csv(secondary_param, "outputs/parameters/secondary_model_mumax_strain_averaged.csv", row.names = FALSE)

# End of saving data into R Project folder 









