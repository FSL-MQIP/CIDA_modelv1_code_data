## ---------------------------Title---------------------------------------------
# CIDA Shelf Life Model v.1.0, sensitivity analysis

## ---------------------------Description---------------------------------------
# Project: CIDA Spinach 

# Description: Sensitivity analysis of the CIDA shelf life model

## ---------------------------Packages------------------------------------------

#   Loading packages
library(dplyr); library(tidyverse); library(ggplot2); library(RColorBrewer); 
library(forcats); library(grDevices); library(extrafont); library(ragg)

## ---------------------------Data----------------------------------------------

#Load in data

b_mumax <- read.csv("sa_outcome/sa_b_mumax_threshold.csv", header = T)
lag_tmin <- read.csv("sa_outcome/sa_lag_tmin_threshold.csv", header = T)
mumax_tmin <- read.csv("sa_outcome/sa_mumax_tmin_threshold.csv", header = T)
mumax <- read.csv("sa_outcome/sa_mumax_threshold.csv", header = T)
n0 <- read.csv("sa_outcome/sa_n0_threshold.csv", header = T)
nmax <- read.csv("sa_outcome/sa_nmax_threshold.csv", header = T)
strain_no <- read.csv("sa_outcome/sa_strain_no_threshold.csv", header = T)
strain_type <- read.csv("sa_outcome/sa_strain_type_threshold.csv", header = T)
lag <- read.csv("sa_outcome/lag_threshold.csv", header = T)

## ---------------------------Data Wrangling------------------------------------

df <- rbind(b_mumax, mumax_tmin, lag_tmin, mumax, n0, nmax, strain_type, lag)

df["baseline"] <- ifelse(df$adj == 1.0, "y", "n")

df$param <- as.factor(df$param)

df <- df %>%
  mutate(percent_fail = prop_spoiled*100) %>%
  group_by(param, day) %>%
  mutate(dev_baseline = percent_fail - percent_fail[baseline == "y"]) %>%
  ungroup()

df <- df %>%
  group_by(day, param) %>%
  mutate(dev_baseline_rel = case_when(
    adj == 0.4 ~ dev_baseline - dev_baseline[adj == 0.6],
    adj == 0.6 ~ dev_baseline - dev_baseline[adj == 0.8],
    adj == 0.8  ~ dev_baseline - dev_baseline[adj == 1.0],
    adj == 1.2 ~ dev_baseline - dev_baseline[adj == 1.0],
    adj == 1.4 ~ dev_baseline - dev_baseline[adj == 1.2],
    adj == 1.6 ~ dev_baseline - dev_baseline[adj == 1.4]
  ))

df$adj <- as.factor(df$adj)

df$adj <- recode(df$adj, "0.4" = "l_60",
                 "0.6" = "l_40",
                 "0.8" = "l_20", 
                 "1" = "b",
                 "1.2" = "u_20", 
                 "1.4" = "u_40",
                 "1.6" = "u_60")

day_18 <- filter(df, day == "18")  %>%
  filter(adj != "b") %>%
  arrange(desc(adj))

day_18_summary <- day_18 %>%
  group_by(param) %>%
  summarise(range = abs(percent_fail[adj == "u_60"] - percent_fail[adj == "l_60"]))

day_18$param <- factor(day_18$param, levels = c("lag_tmin", "fastest_strain_s12_0141", "lag","highest_nmax_s12_0166", "mumax", "mumax_tmin","b_mumax", "n0", "nmax"))

day_18$adj <- factor(day_18$adj, levels = c("l_60", "l_40", "l_20", "u_60", "u_40", "u_20"))

strain_no <- strain_no %>%
  mutate(percent_spoiled = prop_spoiled*100)

strain_no_overall_plotting_data <- strain_no 

strain_no_overall_plotting_data$day <- as.factor(strain_no_overall_plotting_data$day)

## ---------------------------Plotting------------------------------------------

cols <- brewer.pal(6, "PRGn")

day18_tornado_plot <- ggplot(data = day_18, aes(x = param, y = dev_baseline_rel)) + 
  geom_bar(aes(fill = adj), stat = "identity", color = "black") +
  labs(x = "Parameter", y = "Change in packages over the quality threshold, percentage points") + 
  scale_fill_manual(values = cols, labels = c("Decrease, 60%", "Decrease, 40%", "Decrease, 20%", "Increase, 60%", "Increase, 40%", "Increase, 20%")) + 
  guides(fill = guide_legend(title = "Percent change")) + theme_bw() + 
  scale_x_discrete(labels = c(expression("T"["min, "]*lambda), "Strain type:\nFSL S12-0141", expression(lambda), "Strain type:\nFSL S12-0166", 
                              expression(mu["max"]), expression("T"["min, "]*mu[max]), expression("b"[mu[max]]), expression("N"[0]), 
                              expression("N"["max"]))) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2))) 
day18_tornado_plot + coord_flip()

strain_no_overall_plot_spoil <- ggplot(data = strain_no_overall_plotting_data, aes(x = day, y = percent_spoiled, group = adj)) + 
  geom_point(aes(shape = adj)) + scale_shape_manual(values = c(0:6)) + 
  labs(x = "Day", y = "Percentage of packages over the quality threshold") + 
  guides(shape = guide_legend(title = "Number of strains")) + theme_bw() 
strain_no_overall_plot_spoil

## -----------------------------Export------------------------------------------

agg_tiff("tornado_plot_d18_threshold.tiff", width = 150, height = 150, units = "mm", res = 800)
day18_tornado_plot + coord_flip()
dev.off()


agg_tiff("strain_number.tiff", width = 150, height = 150, units = "mm", res = 800)

strain_no_overall_plot_spoil

dev.off()



