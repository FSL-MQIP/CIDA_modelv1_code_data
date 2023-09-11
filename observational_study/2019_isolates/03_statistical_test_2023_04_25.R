##--------------------------------Title-----------------------------------------

#Conducting a statistical test of the association of isolates by the day of isolation

##--------------------------------Description-----------------------------------

#Starting with a chi-sqaured test. If the assumptions for this test are not met, will use a Fisher's exact test

##--------------------------------Packages--------------------------------------

library(dplyr); library(tidyverse); library(vegan)

##--------------------------------Data------------------------------------------

data <- read.csv("data/wrangled/2023_07_09_0919a_isolates_fsl_previd.csv", header = TRUE)

data <- data %>%
  separate(Previous.ID, c("sampling", "sample_name", "day", "media", "isolate_number"), sep = "-") %>%
  separate(subject.titles, c("genus", "species", NA, NA, NA), sep = " ")

at_data <- read.csv("data/raw/16S_AT_assignment.csv", header = TRUE)

##--------------------------------Wrangling-------------------------------------

#Splitting the "Original_ID" column into: sampling, sample replicate, day, test, isolate number

at_data <- separate(at_data, 
                    Original_ID, 
                    into = c("sampling", "sample_replicate", "day", "test", "isolate_number"), 
                    sep = "-")

at_data <- at_data %>%
  filter(Sequenced.isolates %in% data$qseqid)

#Creating a matrix for assessing alpha-diversity, based on 16S AT 
at_table <- as.matrix(table(at_data$day, at_data$AT))

#Categorizing rare genera, i.e., genera which constitute less than 3% of the total isolates (n < 2.88. Rounded up to n < 3)

genera_summary <- data.frame(table(data$genus))

other_genera <- genera_summary %>% 
  filter(Freq < 3) %>% 
  select(Var1)

other_genera <- as.character(other_genera$Var1)

other_genera <- print(other_genera, quote = FALSE)

data$high_freq_genera <- ifelse(data$genus %in% other_genera, "Other", data$genus)

at_freq_table <- at_data %>%
  group_by(AT) %>%
  summarize(count = n())

##--------------------------------Assessing species diversity-------------------

shannon <- diversity(at_table, index = "shannon")

##--------------------------------Test------------------------------------------

test_data <- select(data, c("day", "high_freq_genera"))
test_data$day <- as.factor(test_data$day)
test_data$high_freq_genera <- as.factor(test_data$high_freq_genera)

chisq.test(test_data$day, test_data$high_freq_genera, simulate.p.value = TRUE)$expected
#Expected value is less than 5 in >20% of the cells. So, a Fisher's exact test will be used

test <- fisher.test(table(test_data$day, test_data$high_freq_genera))
##The input table was larger than 2 x 2
##Hybrid = False (default setting). The exact test could be conducted on the full 
# table, so this option was left as false. If hybrid was true, a chi-square test 
# would be conducted if Cochran's conditions were met (these conditions were not met, 
# as found out earlier.)
##Simulate p-value= False. The exact test could be conducted on the full table,
# so there is no need to simulate the p-value
test


