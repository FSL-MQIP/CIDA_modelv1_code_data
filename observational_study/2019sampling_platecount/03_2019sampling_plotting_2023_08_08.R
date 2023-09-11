## --------------------------------Title----------------------------------------
#   Summarizing and plotting the plating data collected during the 2019 sampling

## --------------------------------Description----------------------------------
#   Project: CIDA Spinach 

#  Script description: Plotting the change in APC, GN PC, Y and M over processing and shelf life for the 2019 sampling 

## ----------------------Packages-----------------------------------------------
library(ggplot2); library(dplyr)

## --------------------------------Data-----------------------------------------
data <- read.csv("data/wrangled/0919a_microbial_counts_2023_01_27.csv", header = TRUE)

## --------------------------------Wrangling------------------------------------

data_log <- data %>%
  mutate(log_conc = log10(average_conc)) %>%
  dplyr::select(-average_conc)

#Separating the inprocess and shelf life data 

inprocess_samples <- c("H", "T", "V", "S", "W", "D")

inprocess <- data_log %>%
  filter(sample_id %in% inprocess_samples)

shelflife <- data_log %>%
  filter(!(sample_id %in% inprocess_samples)) %>%
  mutate(day = substr(sample_id, start = 2, stop = nchar(sample_id)))

#Summarizing the average microbial count at different stages of the supply chain, or days of shelf life

summary_inprocess <- inprocess %>%
  group_by(test, sample_id) %>%
  summarize(average = round(mean(log_conc), 2), sd = round(sd(log_conc), 2))


summary_inprocess$sample_id <- factor(summary_inprocess$sample_id, 
                                      levels = c("H", "T", "V", "S", "W", "D"))

summary_inprocess$test <- gsub("Y", "YC", summary_inprocess$test)
summary_inprocess$test <- gsub("M", "MC", summary_inprocess$test)

summary_shelflife <- shelflife %>%
  group_by(test, sample_id, day) %>%
  summarize(average = round(mean(log_conc), 2), sd = round(sd(log_conc), 2))

summary_shelflife$day <- factor(summary_shelflife$day, 
                                levels = c("0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20"))

colnames(summary_shelflife)[1] <- "Test"

summary_shelflife$Test <- gsub("Y", "YC", summary_shelflife$Test)
summary_shelflife$Test <- gsub("M", "MC", summary_shelflife$Test)

## --------------------------------Plotting-----------------------------------------

date <- Sys.Date()
date <- gsub("-", "_", date)

#In process
png(paste0("outputs/", date, "_microbial_counts_over_processing.png"), width = 150, height = 150, units = "mm", res = 1200)
ggplot(data = summary_inprocess, mapping = aes(x=sample_id, y=average, group = test)) + 
  labs(x = "Step in the Supply Chain", y = expression("Microbial concentration, log"[10]*"CFU/g")) + 
  scale_y_continuous(limits = c(1.5, 7.5)) + geom_point(position = position_dodge(0.3), aes(shape = test, group = test) ) + 
  geom_errorbar(aes(ymin = average - 0.5*sd, ymax = average + 0.5*sd), width = 0.01, position = position_dodge(0.3)) +
  scale_shape_manual(values = c(0, 1, 2, 5, 6)) + labs(shape = "Test") + theme_bw()
dev.off()

#Shelf life
png(paste0("outputs/", date, "_microbial_counts_over_shelf_life.png"), width = 150, height = 150, units = "mm", res = 1200)
ggplot(data = summary_shelflife, mapping = aes(x = day, y = average, group = Test)) +  
  geom_line(data = summary_shelflife, aes(x = day, y = average, group = Test, linetype = Test)) + 
  geom_point(position = position_dodge(0.3), aes(shape = Test)) + 
  geom_errorbar(aes(ymin = average - 0.5*sd, ymax = average + 0.5*sd), width = 0.01, position = position_dodge(0.3)) +
  labs(x = "Day", y = expression("Microbial concentration, log"[10]*"CFU/g")) +
  coord_fixed() + theme_bw() + 
  scale_shape_manual(values = c(0, 1, 2, 5, 6)) 
dev.off()

