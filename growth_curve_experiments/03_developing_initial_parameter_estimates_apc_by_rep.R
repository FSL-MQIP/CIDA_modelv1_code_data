##  ----------------------------Title-------------------------------------------
#   Getting initial parameter estimates for SM-APC

##  --------------------------Description---------------------------------------
#   Project: CIDA Spinach 

#   Script description: Initial parameter estimates, for primary growth models, are assessed for all SM-APC replicates 

##  --------------------------Packages------------------------------------------
#   Loading packages
library(tidyverse); library(dplyr); library(minpack.lm)

##  --------------------------Data----------------------------------------------
#Read in data 
apc_strain_data <- read.csv("data/apc_strain_averaged_wrangled_02.csv", header = TRUE)

## ---------------------------Data Wrangling------------------------------------
#Data is split into strain growth data and APC growth data 
strain_data <- apc_strain_data %>%
  dplyr::filter(media == "BHI")

#Subset strain data by isolate and rep
b116_6_1 <- filter(strain_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B2")
b116_6_2 <- filter(strain_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B3")
b116_6_3 <- filter(strain_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "6" & batch == "B6")

b132_6_1 <- filter(strain_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B2")
b132_6_2 <- filter(strain_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B3")
b132_6_3 <- filter(strain_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "6" & batch == "B6")

b141_6_1 <- filter(strain_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B2")
b141_6_2 <- filter(strain_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B3")
b141_6_3 <- filter(strain_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "6" & batch == "B6")

b166_6_1 <- filter(strain_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B4")
b166_6_2 <- filter(strain_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B5")
b166_6_3 <- filter(strain_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "6" & batch == "B7")

b180_6_1 <- filter(strain_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B4")
b180_6_2 <- filter(strain_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B5")
b180_6_3 <- filter(strain_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "6" & batch == "B7")

b184_6_1 <- filter(strain_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B4")
b184_6_2 <- filter(strain_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B5")
b184_6_3 <- filter(strain_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "6" & batch == "B7")

b116_10 <- filter(strain_data, isolate ==  c("S12-0116")) %>%
  filter(temperature == "10")
b132_10 <- filter(strain_data, isolate ==  c("S12-0132")) %>%
  filter(temperature == "10")
b141_10 <- filter(strain_data, isolate ==  c("S12-0141")) %>%
  filter(temperature == "10")
b166_10 <- filter(strain_data, isolate ==  c("S12-0166")) %>%
  filter(temperature == "10")
b180_10 <- filter(strain_data, isolate ==  c("S12-0180")) %>%
  filter(temperature == "10")
b184_10 <- filter(strain_data, isolate ==  c("S12-0184")) %>%
  filter(temperature == "10")

#Making a dataframe to store initial parameter estimates 
no_of_isolates <- 6
no_of_temps <- 2
no_of_reps <- 3
isolates_names <- factor(c("S12-0116", "S12-0132", "S12-0141", "S12-0166", "S12-0180", "S12-0184"))
param_estimates <- data.frame(matrix(nrow = no_of_isolates*4, ncol = 7))
colnames(param_estimates) <- c("isolate", "n0", "lag", "mumax", "nmax", "temp", "rep")
param_estimates$isolate <- c(as.character(rep(isolates_names, each = 3)), as.character(isolates_names))
param_estimates$temp <- c(rep("6", times = 18), rep("10", times = 6))
param_estimates$rep <- c(rep(1:3, times = 6), rep(1, times = 6))

#End of data wrangling

##----------------------------S12-0116, 6C--------------------------------------
#S12-0116, rep 1
#n0 
param_estimates[1, 2] <- b116_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave)) 

#lag 
for(i in 1:nrow(b116_6_1)){
  if((b116_6_1$log_average_wrangled_conc[i] - b116_6_1$log_average_wrangled_conc[1])/b116_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[1,3] <- b116_6_1$day[i]
    break
  }
}

#mu 
for(i in 1:nrow(b116_6_1)){
  if((b116_6_1$log_average_wrangled_conc[i] - b116_6_1$log_average_wrangled_conc[1])/b116_6_1$log_average_wrangled_conc[1] > 0.15){
    b116_6_1_x1 <- b116_6_1$day[i]
    b116_6_1_y1 <- b116_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b116_6_1)){
  if(round(b116_6_1$log_average_wrangled_conc[i]) == round(max(b116_6_1$log_average_wrangled_conc))){
    b116_6_1_x2 <- b116_6_1$day[i]
    b116_6_1_y2 <- b116_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[1,4]  <- (b116_6_1_y2 - b116_6_1_y1)/(b116_6_1_x2 - b116_6_1_x1)

#nmax 
param_estimates[1,5] <- b116_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0116, rep 2
#n0 
param_estimates[2, 2] <- b116_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave)) 

#lag 
for(i in 1:nrow(b116_6_2)){
  if((b116_6_2$log_average_wrangled_conc[i] - b116_6_2$log_average_wrangled_conc[1])/b116_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[2,3] <- b116_6_2$day[i]
    break
  }
}

#mu 
for(i in 1:nrow(b116_6_2)){
  if((b116_6_2$log_average_wrangled_conc[i] - b116_6_2$log_average_wrangled_conc[1])/b116_6_2$log_average_wrangled_conc[1] > 0.15){
    b116_6_2_x1 <- b116_6_2$day[i]
    b116_6_2_y1 <- b116_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b116_6_2)){
  if(round(b116_6_2$log_average_wrangled_conc[i]) == round(max(b116_6_2$log_average_wrangled_conc))){
    b116_6_2_x2 <- b116_6_2$day[i]
    b116_6_2_y2 <- b116_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[2,4]  <- (b116_6_2_y2 - b116_6_2_y1)/(b116_6_2_x2 - b116_6_2_x1)

#nmax 
param_estimates[2,5] <- b116_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0116, rep 3
#n0 
param_estimates[3, 2] <- b116_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave)) 

#lag
for(i in 1:nrow(b116_6_3)){
  if((b116_6_3$log_average_wrangled_conc[i] - b116_6_3$log_average_wrangled_conc[1])/b116_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[3,3] <- b116_6_3$day[i]
    break
  }
}

#mu 
for(i in 1:nrow(b116_6_3)){
  if((b116_6_3$log_average_wrangled_conc[i] - b116_6_3$log_average_wrangled_conc[1])/b116_6_3$log_average_wrangled_conc[1] > 0.15){
    b116_6_3_x1 <- b116_6_3$day[i]
    b116_6_3_y1 <- b116_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b116_6_3)){
  if(round(b116_6_3$log_average_wrangled_conc[i]) == round(max(b116_6_3$log_average_wrangled_conc))){
    b116_6_3_x2 <- b116_6_3$day[i]
    b116_6_3_y2 <- b116_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[3,4]  <- (b116_6_3_y2 - b116_6_3_y1)/(b116_6_3_x2 - b116_6_3_x1)

#nmax 
param_estimates[3,5] <- b116_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 


##----------------------------S12-0132, 6C--------------------------------------
#S12-0132, rep 1
#n0 
param_estimates[4,2] <- b132_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b132_6_1)){
  if((b132_6_1$log_average_wrangled_conc[i] - b132_6_1$log_average_wrangled_conc[1])/b132_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[4,3] <- print(b132_6_1$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b132_6_1)){
  if((b132_6_1$log_average_wrangled_conc[i] - b132_6_1$log_average_wrangled_conc[1])/b132_6_1$log_average_wrangled_conc[1] > 0.15){
    b132_6_1_x1 <- b132_6_1$day[i]
    b132_6_1_y1 <- b132_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b132_6_1)){
  if(round(b132_6_1$log_average_wrangled_conc[i]) == round(max(b132_6_1$log_average_wrangled_conc))){
    b132_6_1_x2 <- b132_6_1$day[i]
    b132_6_1_y2 <- b132_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[4,4]  <- (b132_6_1_y2 - b132_6_1_y1)/(b132_6_1_x2 - b132_6_1_x1)

#nmax 
param_estimates[4,5] <- b132_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0132, rep 2
#n0 
param_estimates[5,2] <- b132_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b132_6_2)){
  if((b132_6_2$log_average_wrangled_conc[i] - b132_6_2$log_average_wrangled_conc[1])/b132_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[5,3] <- print(b132_6_2$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b132_6_2)){
  if((b132_6_2$log_average_wrangled_conc[i] - b132_6_2$log_average_wrangled_conc[1])/b132_6_2$log_average_wrangled_conc[1] > 0.15){
    b132_6_2_x1 <- b132_6_2$day[i]
    b132_6_2_y1 <- b132_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b132_6_2)){
  if(round(b132_6_2$log_average_wrangled_conc[i]) == round(max(b132_6_2$log_average_wrangled_conc))){
    b132_6_2_x2 <- b132_6_2$day[i]
    b132_6_2_y2 <- b132_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[5,4]  <- (b132_6_2_y2 - b132_6_2_y1)/(b132_6_2_x2 - b132_6_2_x1)

#nmax 
param_estimates[5,5] <- b132_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0132, rep 3
#n0 
param_estimates[6,2] <- b132_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b132_6_3)){
  if((b132_6_3$log_average_wrangled_conc[i] - b132_6_3$log_average_wrangled_conc[1])/b132_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[6,3] <- print(b132_6_3$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b132_6_3)){
  if((b132_6_3$log_average_wrangled_conc[i] - b132_6_3$log_average_wrangled_conc[1])/b132_6_3$log_average_wrangled_conc[1] > 0.15){
    b132_6_3_x1 <- b132_6_3$day[i]
    b132_6_3_y1 <- b132_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b132_6_3)){
  if(round(b132_6_3$log_average_wrangled_conc[i]) == round(max(b132_6_3$log_average_wrangled_conc))){
    b132_6_3_x2 <- b132_6_3$day[i]
    b132_6_3_y2 <- b132_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[6,4]  <- (b132_6_3_y2 - b132_6_3_y1)/(b132_6_3_x2 - b132_6_3_x1)

#nmax 
param_estimates[6,5] <- b132_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0141, 6C--------------------------------------

#S12-0141, rep 1
#n0 
param_estimates[7,2] <- b141_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b141_6_1)){
  if((b141_6_1$log_average_wrangled_conc[i] - b141_6_1$log_average_wrangled_conc[1])/b141_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[7,3] <- print(b141_6_1$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b141_6_1)){
  if((b141_6_1$log_average_wrangled_conc[i] - b141_6_1$log_average_wrangled_conc[1])/b141_6_1$log_average_wrangled_conc[1] > 0.15){
    b141_6_1_x1 <- b141_6_1$day[i]
    b141_6_1_y1 <- b141_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b141_6_1)){
  if(round(b141_6_1$log_average_wrangled_conc[i]) == round(max(b141_6_1$log_average_wrangled_conc))){
    b141_6_1_x2 <- b141_6_1$day[i]
    b141_6_1_y2 <- b141_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[7,4]  <- (b141_6_1_y2 - b141_6_1_y1)/(b141_6_1_x2 - b141_6_1_x1)

#nmax 
param_estimates[7,5] <- b141_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0141, rep 2
#n0 
param_estimates[8,2] <- b141_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b141_6_2)){
  if((b141_6_2$log_average_wrangled_conc[i] - b141_6_2$log_average_wrangled_conc[1])/b141_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[8,3] <- print(b141_6_2$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b141_6_2)){
  if((b141_6_2$log_average_wrangled_conc[i] - b141_6_2$log_average_wrangled_conc[1])/b141_6_2$log_average_wrangled_conc[1] > 0.15){
    b141_6_2_x1 <- b141_6_2$day[i]
    b141_6_2_y1 <- b141_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b141_6_2)){
  if(round(b141_6_2$log_average_wrangled_conc[i]) == round(max(b141_6_2$log_average_wrangled_conc))){
    b141_6_2_x2 <- b141_6_2$day[i]
    b141_6_2_y2 <- b141_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[8,4]  <- (b141_6_2_y2 - b141_6_2_y1)/(b141_6_2_x2 - b141_6_2_x1)

#nmax 
param_estimates[8,5] <- b141_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0141, rep 3
#n0 
param_estimates[9,2] <- b141_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b141_6_3)){
  if((b141_6_3$log_average_wrangled_conc[i] - b141_6_3$log_average_wrangled_conc[1])/b141_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[9,3] <- print(b141_6_3$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b141_6_3)){
  if((b141_6_3$log_average_wrangled_conc[i] - b141_6_3$log_average_wrangled_conc[1])/b141_6_3$log_average_wrangled_conc[1] > 0.15){
    b141_6_3_x1 <- b141_6_3$day[i]
    b141_6_3_y1 <- b141_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b141_6_3)){
  if(round(b141_6_3$log_average_wrangled_conc[i]) == round(max(b141_6_3$log_average_wrangled_conc))){
    b141_6_3_x2 <- b141_6_3$day[i]
    b141_6_3_y2 <- b141_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[9,4]  <- (b141_6_3_y2 - b141_6_3_y1)/(b141_6_3_x2 - b141_6_3_x1)

#nmax 
param_estimates[9,5] <- b141_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0166, 6C--------------------------------------
#S12-0166, rep 1
#Initial parameter estimates for S12-0166, 6C
#n0 
param_estimates[10,2] <- b166_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b166_6_1)){
  if((b166_6_1$log_average_wrangled_conc[i] - b166_6_1$log_average_wrangled_conc[1])/b166_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[10,3] <- print(b166_6_1$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b166_6_1)){
  if((b166_6_1$log_average_wrangled_conc[i] - b166_6_1$log_average_wrangled_conc[1])/b166_6_1$log_average_wrangled_conc[1] > 0.15){
    b166_6_1_x1 <- b166_6_1$day[i]
    b166_6_1_y1 <- b166_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b166_6_1)){
  if(round(b166_6_1$log_average_wrangled_conc[i]) == round(max(b166_6_1$log_average_wrangled_conc))){
    b166_6_1_x2 <- b166_6_1$day[i]
    b166_6_1_y2 <- b166_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[10,4]  <- (b166_6_1_y2 - b166_6_1_y1)/(b166_6_1_x2 - b166_6_1_x1)

#nmax 
param_estimates[10,5] <- b166_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0166, rep 2
#Initial parameter estimates for S12-0166, 6C
#n0 
param_estimates[11,2] <- b166_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b166_6_2)){
  if((b166_6_2$log_average_wrangled_conc[i] - b166_6_2$log_average_wrangled_conc[1])/b166_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[11,3] <- print(b166_6_2$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b166_6_2)){
  if((b166_6_2$log_average_wrangled_conc[i] - b166_6_2$log_average_wrangled_conc[1])/b166_6_2$log_average_wrangled_conc[1] > 0.15){
    b166_6_2_x1 <- b166_6_2$day[i]
    b166_6_2_y1 <- b166_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b166_6_2)){
  if(round(b166_6_2$log_average_wrangled_conc[i]) == round(max(b166_6_2$log_average_wrangled_conc))){
    b166_6_2_x2 <- b166_6_2$day[i]
    b166_6_2_y2 <- b166_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[11,4]  <- (b166_6_2_y2 - b166_6_2_y1)/(b166_6_2_x2 - b166_6_2_x1)

#nmax 
param_estimates[11,5] <- b166_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0166, rep 3
#Initial parameter estimates for S12-0166, 6C
#n0 
param_estimates[12,2] <- b166_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b166_6_3)){
  if((b166_6_3$log_average_wrangled_conc[i] - b166_6_3$log_average_wrangled_conc[1])/b166_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[12,3] <- print(b166_6_3$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b166_6_3)){
  if((b166_6_3$log_average_wrangled_conc[i] - b166_6_3$log_average_wrangled_conc[1])/b166_6_3$log_average_wrangled_conc[1] > 0.15){
    b166_6_3_x1 <- b166_6_3$day[i]
    b166_6_3_y1 <- b166_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b166_6_3)){
  if(round(b166_6_3$log_average_wrangled_conc[i]) == round(max(b166_6_3$log_average_wrangled_conc))){
    b166_6_3_x2 <- b166_6_3$day[i]
    b166_6_3_y2 <- b166_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[12,4]  <- (b166_6_3_y2 - b166_6_3_y1)/(b166_6_3_x2 - b166_6_3_x1)

#nmax 
param_estimates[12,5] <- b166_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 


##----------------------------S12-0180, 6C--------------------------------------

#S12-0180, rep 1
#n0 
param_estimates[13,2] <- b180_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b180_6_1)){
  if((b180_6_1$log_average_wrangled_conc[i] - b180_6_1$log_average_wrangled_conc[1])/b180_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[13,3] <- print(b180_6_1$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b180_6_1)){
  if((b180_6_1$log_average_wrangled_conc[i] - b180_6_1$log_average_wrangled_conc[1])/b180_6_1$log_average_wrangled_conc[1] > 0.15){
    b180_6_1_x1 <- b180_6_1$day[i]
    b180_6_1_y1 <- b180_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b180_6_1)){
  if(round(b180_6_1$log_average_wrangled_conc[i]) == round(max(b180_6_1$log_average_wrangled_conc))){
    b180_6_1_x2 <- b180_6_1$day[i]
    b180_6_1_y2 <- b180_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[13,4]  <- (b180_6_1_y2 - b180_6_1_y1)/(b180_6_1_x2 - b180_6_1_x1)

#nmax
param_estimates[13,5] <- b180_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0180, rep 2
#n0 
param_estimates[14,2] <- b180_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag
for(i in 1:nrow(b180_6_2)){
  if((b180_6_2$log_average_wrangled_conc[i] - b180_6_2$log_average_wrangled_conc[1])/b180_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[14,3] <- print(b180_6_2$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b180_6_2)){
  if((b180_6_2$log_average_wrangled_conc[i] - b180_6_2$log_average_wrangled_conc[1])/b180_6_2$log_average_wrangled_conc[1] > 0.15){
    b180_6_2_x1 <- b180_6_2$day[i]
    b180_6_2_y1 <- b180_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b180_6_2)){
  if(round(b180_6_2$log_average_wrangled_conc[i]) == round(max(b180_6_2$log_average_wrangled_conc))){
    b180_6_2_x2 <- b180_6_2$day[i]
    b180_6_2_y2 <- b180_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[14,4]  <- (b180_6_2_y2 - b180_6_2_y1)/(b180_6_2_x2 - b180_6_2_x1)

#nmax
param_estimates[14,5] <- b180_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave))

#S12-0180, rep 3
#n0 
param_estimates[15,2] <- b180_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b180_6_3)){
  if((b180_6_3$log_average_wrangled_conc[i] - b180_6_3$log_average_wrangled_conc[1])/b180_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[15,3] <- print(b180_6_3$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b180_6_3)){
  if((b180_6_3$log_average_wrangled_conc[i] - b180_6_3$log_average_wrangled_conc[1])/b180_6_3$log_average_wrangled_conc[1] > 0.15){
    b180_6_3_x1 <- b180_6_3$day[i]
    b180_6_3_y1 <- b180_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b180_6_3)){
  if(round(b180_6_3$log_average_wrangled_conc[i]) == round(max(b180_6_3$log_average_wrangled_conc))){
    b180_6_3_x2 <- b180_6_3$day[i]
    b180_6_3_y2 <- b180_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[15,4]  <- (b180_6_3_y2 - b180_6_3_y1)/(b180_6_3_x2 - b180_6_3_x1)

#nmax
param_estimates[15,5] <- b180_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave))

##----------------------------S12-0184, 6C--------------------------------------

#S12-0184, rep 1
#n0 
param_estimates[16,2] <- b184_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag
for(i in 1:nrow(b184_6_1)){
  if((b184_6_1$log_average_wrangled_conc[i] - b184_6_1$log_average_wrangled_conc[1])/b184_6_1$log_average_wrangled_conc[1] > 0.15){
    param_estimates[16,3] <- print(b184_6_1$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b184_6_1)){
  if((b184_6_1$log_average_wrangled_conc[i] - b184_6_1$log_average_wrangled_conc[1])/b184_6_1$log_average_wrangled_conc[1] > 0.15){
    b184_6_1_x1 <- b184_6_1$day[i]
    b184_6_1_y1 <- b184_6_1$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b184_6_1)){
  if(round(b184_6_1$log_average_wrangled_conc[i]) == round(max(b184_6_1$log_average_wrangled_conc))){
    b184_6_1_x2 <- b184_6_1$day[i]
    b184_6_1_y2 <- b184_6_1$log_average_wrangled_conc[i]
    break
  }}

param_estimates[16,4]  <- (b184_6_1_y2 - b184_6_1_y1)/(b184_6_1_x2 - b184_6_1_x1)

#nmax
param_estimates[16,5] <- b184_6_1 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0184, rep 2
#n0 
param_estimates[17,2] <- b184_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b184_6_2)){
  if((b184_6_2$log_average_wrangled_conc[i] - b184_6_2$log_average_wrangled_conc[1])/b184_6_2$log_average_wrangled_conc[1] > 0.15){
    param_estimates[17,3] <- print(b184_6_2$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b184_6_2)){
  if((b184_6_2$log_average_wrangled_conc[i] - b184_6_2$log_average_wrangled_conc[1])/b184_6_2$log_average_wrangled_conc[1] > 0.15){
    b184_6_2_x1 <- b184_6_2$day[i]
    b184_6_2_y1 <- b184_6_2$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b184_6_2)){
  if(round(b184_6_2$log_average_wrangled_conc[i]) == round(max(b184_6_2$log_average_wrangled_conc))){
    b184_6_2_x2 <- b184_6_2$day[i]
    b184_6_2_y2 <- b184_6_2$log_average_wrangled_conc[i]
    break
  }}

param_estimates[17,4]  <- (b184_6_2_y2 - b184_6_2_y1)/(b184_6_2_x2 - b184_6_2_x1)

#nmax
param_estimates[17,5] <- b184_6_2 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#S12-0184, rep 3
#n0 
param_estimates[18,2] <- b184_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag
for(i in 1:nrow(b184_6_3)){
  if((b184_6_3$log_average_wrangled_conc[i] - b184_6_3$log_average_wrangled_conc[1])/b184_6_3$log_average_wrangled_conc[1] > 0.15){
    param_estimates[18,3] <- print(b184_6_3$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b184_6_3)){
  if((b184_6_3$log_average_wrangled_conc[i] - b184_6_3$log_average_wrangled_conc[1])/b184_6_3$log_average_wrangled_conc[1] > 0.15){
    b184_6_3_x1 <- b184_6_3$day[i]
    b184_6_3_y1 <- b184_6_3$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b184_6_3)){
  if(round(b184_6_3$log_average_wrangled_conc[i]) == round(max(b184_6_3$log_average_wrangled_conc))){
    b184_6_3_x2 <- b184_6_3$day[i]
    b184_6_3_y2 <- b184_6_3$log_average_wrangled_conc[i]
    break
  }}

param_estimates[18,4]  <- (b184_6_3_y2 - b184_6_3_y1)/(b184_6_3_x2 - b184_6_3_x1)

#nmax
param_estimates[18,5] <- b184_6_3 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0116, 10C-------------------------------------

#Initial parameter estimates for S12-0116, 10C 
#n0 
param_estimates[19, 2] <- b116_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave)) 

#lag 
for(i in 1:nrow(b116_10)){
  if((b116_10$log_average_wrangled_conc[i] - b116_10$log_average_wrangled_conc[1])/b116_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[19,3] <- print(b116_10$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b116_10)){
  if((b116_10$log_average_wrangled_conc[i] - b116_10$log_average_wrangled_conc[1])/b116_10$log_average_wrangled_conc[1] > 0.15){
    b116_10_x1 <- b116_10$day[i]
    b116_10_y1 <- b116_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b116_10)){
  if(round(b116_10$log_average_wrangled_conc[i]) == round(max(b116_10$log_average_wrangled_conc))){
    b116_10_x2 <- b116_10$day[i]
    b116_10_y2 <- b116_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[19,4]  <- (b116_10_y2 - b116_10_y1)/(b116_10_x2 - b116_10_x1)

#nmax 
param_estimates[19,5] <- b116_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0132, 10C-------------------------------------
#n0 
param_estimates[20,2] <- b132_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag
for(i in 1:nrow(b132_10)){
  if((b132_10$log_average_wrangled_conc[i] - b132_10$log_average_wrangled_conc[1])/b132_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[20,3] <- print(b132_10$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b132_10)){
  if((b132_10$log_average_wrangled_conc[i] - b132_10$log_average_wrangled_conc[1])/b132_10$log_average_wrangled_conc[1] > 0.15){
    b132_10_x1 <- b132_10$day[i]
    b132_10_y1 <- b132_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b132_10)){
  if(round(b132_10$log_average_wrangled_conc[i]) == round(max(b132_10$log_average_wrangled_conc))){
    b132_10_x2 <- b132_10$day[i]
    b132_10_y2 <- b132_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[20,4]  <- (b132_10_y2 - b132_10_y1)/(b132_10_x2 - b132_10_x1)

#nmax 
param_estimates[20,5] <- b132_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0141, 10C-------------------------------------
#n0 
param_estimates[21,2] <- b141_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b141_10)){
  if((b141_10$log_average_wrangled_conc[i] - b141_10$log_average_wrangled_conc[1])/b141_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[21,3] <- print(b141_10$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b141_10)){
  if((b141_10$log_average_wrangled_conc[i] - b141_10$log_average_wrangled_conc[1])/b141_10$log_average_wrangled_conc[1] > 0.15){
    b141_10_x1 <- b141_10$day[i]
    b141_10_y1 <- b141_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b141_10)){
  if(round(b141_10$log_average_wrangled_conc[i]) == round(max(b141_10$log_average_wrangled_conc))){
    b141_10_x2 <- b141_10$day[i]
    b141_10_y2 <- b141_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[21,4]  <- (b141_10_y2 - b141_10_y1)/(b141_10_x2 - b141_10_x1)

#nmax 
param_estimates[21,5] <- b141_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0166, 10C-------------------------------------
#n0 
param_estimates[22,2] <- b166_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag (changed to 1, based on visual assessment of a figure of counts by day for this isolate/replicate)
for(i in 1:nrow(b166_10)){
  if((b166_10$log_average_wrangled_conc[i] - b166_10$log_average_wrangled_conc[1])/b166_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[22,3] <- print(b166_10$day[i])
    break
  }
}

param_estimates[22,3] <- "1"

#mu 
for(i in 1:nrow(b166_10)){
  if((b166_10$log_average_wrangled_conc[i] - b166_10$log_average_wrangled_conc[1])/b166_10$log_average_wrangled_conc[1] > 0.15){
    b166_10_x1 <- b166_10$day[i]
    b166_10_y1 <- b166_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b166_10)){
  if(round(b166_10$log_average_wrangled_conc[i]) == round(max(b166_10$log_average_wrangled_conc))){
    b166_10_x2 <- b166_10$day[i]
    b166_10_y2 <- b166_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[22,4]  <- (b166_10_y2 - b166_10_y1)/(b166_10_x2 - b166_10_x1)

#nmax 
param_estimates[22,5] <- b166_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0180, 10C-------------------------------------
#n0 
param_estimates[23,2] <- b180_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag (replaced with 2, based on visual assessment)
for(i in 1:nrow(b180_10)){
  if((b180_10$log_average_wrangled_conc[i] - b180_10$log_average_wrangled_conc[1])/b180_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[23,3] <- print(b180_10$day[i])
    break
  }
}

param_estimates[23,3] <- 2

#mu 
for(i in 1:nrow(b180_10)){
  if((b180_10$log_average_wrangled_conc[i] - b180_10$log_average_wrangled_conc[1])/b180_10$log_average_wrangled_conc[1] > 0.15){
    b180_10_x1 <- b180_10$day[i]
    b180_10_y1 <- b180_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b180_10)){
  if(round(b180_10$log_average_wrangled_conc[i]) == round(max(b180_10$log_average_wrangled_conc))){
    b180_10_x2 <- b180_10$day[i]
    b180_10_y2 <- b180_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[23,4]  <- (b180_10_y2 - b180_10_y1)/(b180_10_x2 - b180_10_x1)

#nmax
param_estimates[23,5] <- b180_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

##----------------------------S12-0184, 10C-------------------------------------
#n0 
param_estimates[24,2] <- b184_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(min(ave))

#lag 
for(i in 1:nrow(b184_10)){
  if((b184_10$log_average_wrangled_conc[i] - b184_10$log_average_wrangled_conc[1])/b184_10$log_average_wrangled_conc[1] > 0.15){
    param_estimates[24,3] <- print(b184_10$day[i])
    break
  }
}

#mu 
for(i in 1:nrow(b184_10)){
  if((b184_10$log_average_wrangled_conc[i] - b184_10$log_average_wrangled_conc[1])/b184_10$log_average_wrangled_conc[1] > 0.15){
    b184_10_x1 <- b184_10$day[i]
    b184_10_y1 <- b184_10$log_average_wrangled_conc[i]
    break
  }}

for(i in 1:nrow(b184_10)){
  if(round(b184_10$log_average_wrangled_conc[i]) == round(max(b184_10$log_average_wrangled_conc))){
    b184_10_x2 <- b184_10$day[i]
    b184_10_y2 <- b184_10$log_average_wrangled_conc[i]
    break
  }}

param_estimates[24,4]  <- (b184_10_y2 - b184_10_y1)/(b184_10_x2 - b184_10_x1)

#nmax
param_estimates[24,5] <- b184_10 %>%
  group_by(day) %>%
  summarise(ave = mean(log_average_wrangled_conc)) %>%
  summarise(max(ave)) 

#End of Initial Parameter Estimates, 10C

##------------------------Exporting Initial Model Parameters--------------------

#Push the data back to the R project 

#write.csv(param_estimates, "outputs/parameters/initial_parameter_estimates_apc_by_rep.csv", row.names = FALSE)

