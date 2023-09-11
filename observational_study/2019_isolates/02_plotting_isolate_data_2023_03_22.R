##--------------------------------Title-----------------------------------------

#Plotting isolate sequencing results

##--------------------------------Description-----------------------------------

#Preparing a barplot of genus-level identity by day of isolation. A color and
#grayscale version of the plot will be prepared

##--------------------------------Packages--------------------------------------

library(tidyverse); library(ggplot2); library(stringi);

##--------------------------------Data------------------------------------------

data <- read.csv("data/wrangled/2023_07_09_0919a_isolates_fsl_previd.csv", header = TRUE)

##--------------------------------Wrangling-------------------------------------

data <- data %>%
  separate(Previous.ID, c("sampling", "sample_name", "day", "media", "isolate_number"), sep = "-") %>%
  separate(subject.titles, c("genus", "species", NA, NA, NA, NA), sep = " ", extra = "drop")

data$day <- stri_replace_all_regex(data$day, 
                                   pattern = c("D0", "D4", "D10", "D20"),
                                   replacement = c(0, 4, 10, 20),
                                   vectorize_all = FALSE)

data$day <- factor(data$day, levels = c("0", "4", "10", "20"))

freq <- data.frame(table(data$genus, data$day))

colnames(freq) <- c("genus", "day", "freq")

freq_plot_data <- freq %>%
  filter(freq != 0)

freq_table <- table(data$genus, data$day)
##--------------------------------Plotting--------------------------------------

#Color plot
ggplot(data = freq_plot_data, mapping = aes(x = day, y = freq, fill = genus, label = freq)) + 
  geom_col(position = "stack", color = "black") +
  theme_classic() + ylab("Number of isolates") + xlab("Day") +
  labs(fill = "Genus") + scale_fill_manual(values = c("#AE76A3","#1965B0", "#5289C7", "#1E9F02", 
                                                      "#CAE0AB", "#F7F056", "#F6C141","#F1932D", "#E8601C", 
                                                      "#AFADAD")) +
  geom_text(position = position_stack(vjust = 0.5))

#Greyscale plot
angle1 <- rep(c(0, 45,45,135), length.out=10)
angle2 <- rep(c(0, 45,135,135), length.out=10)
density1 <- seq(5,20,length.out=10)
density2 <- seq(5,20,length.out=10)
col <- 1 

barplot(freq_table, col = col, angle = angle1, density = density1)
barplot(freq_table, col = col, angle = angle2, density = density2)
legend("topright", legend=1:10, ncol =1, fill=TRUE, col=col, angle=angle1, density=density1)
par(bg="transparent")
legend("topright", legend=1:10, ncol =1, fill=TRUE, col=col, angle=angle2, density=density2)
par(op)


##--------------------------------Export----------------------------------------

date <- Sys.Date()
date <- gsub("-", "_", date)

png(paste0("outputs/", date, "_microbial_population_over_shelf_life.png"), width = 150, height = 150, units = "mm", res = 1200)
ggplot(data = freq_plot_data, mapping = aes(x = day, y = freq, fill = genus, label = freq)) + 
  geom_col(position = "stack", color = "black") +
  theme_classic() + ylab("Number of isolates") + xlab("Day") +
  labs(fill = "Genus") + scale_fill_manual(values = c("#AE76A3","#1965B0", "#5289C7", "#1E9F02", 
                                                      "#CAE0AB", "#F7F056", "#F6C141","#F1932D", "#E8601C", 
                                                      "#AFADAD")) +
  geom_text(position = position_stack(vjust = 0.5))
dev.off()
