##  ----------------------------Title-------------------------------------------
#Assigning ST frequency for the CIDA isolate data from the 2019 sampling
##  --------------------------Description---------------------------------------
#   Project: CIDA Project 

#   Script description: The ST frequency will be used to assign the strains (used in the CIDA growth experiments) to packages of spinach in the Monte Carlo simulation
##  --------------------------Packages------------------------------------------
library(dplyr)
library(reshape2)

##----------------------------Data----------------------------------------------
data <- read.csv("data/raw/distance_matrix.csv", header = T, check.names = F) 

##----------------------------Data Wrangling------------------------------------
colnames(data)[1] <- "isolates"

cida_growth_isolates <- c("S12-0116:0919aM9-0536", "S12-0132:0919aM9-0540", "S12-0141:0919aM9-0606", "S12-0166:0919aS12-0004", "S12-0180:0919aS12-0063", "S12-0184:0919aS12-0074")

data_2 <- melt(data)

data_3 <- data_2 %>%
  filter(!(isolates %in% cida_growth_isolates)) %>%
  filter(variable %in% cida_growth_isolates)

colnames(data_3) <- c("all_isolates", "growth_isolates", "percent_identity")

data_4 <- data_3 %>%
  group_by(all_isolates) %>%
  arrange(desc(percent_identity)) %>%
  slice(1)

data_5 <- data_4 %>%
  group_by(growth_isolates) %>%
  summarise(count = n()) %>%
  mutate(frequency = count/sum(count))

##-----------------------------Export-------------------------------------------

#Push the data back to the R project 
date <- Sys.Date()
date <- gsub("-", "_", date)

write.csv(data_5, "outputs/parameters/st_frequency.csv", row.names = FALSE)


