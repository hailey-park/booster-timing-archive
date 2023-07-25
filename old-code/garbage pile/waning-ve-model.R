###################################################################################################
#Title: Waning VE Model
#Author: Hailey Park
#Date: March 25, 2023
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo/Booster Timing")

#Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)
library(scales)
library(data.table)
library(lme4)
library(merTools)


#Load in dataset
waning_data_age_strat_alldose <- read.csv("data/waning_data_absolute_new_agestrat_alldose.csv") #USING THIS ONE
waning_relative_4dose <- read.csv("data/relative-ve-4dose-lin.csv")
waning_data_95UI <- read.csv("data/waning_data_absolute_new_agestrat_alldose_95UI.csv")

#Clean
waning_relative_4dose <- waning_relative_4dose %>% mutate(months = log(Month)) %>% rename(relative_ve = ve)


#Separate lower and upper 95% UI estimates
lower_95UI <- waning_data_95UI %>% filter(VE.Type == "lower")
upper_95UI <- waning_data_95UI %>% filter(VE.Type == "upper")


write.csv(rbind(read.csv("ve_waning_predictions_adj_upper.csv")[,-1], read.csv("ve_waning_predictions_adj_lower.csv")[,-1]), "ve_waning_predictions_adj_95UI.csv")


#Reformat data to long and clean data
long_data <- melt(waning_data_age_strat_alldose) %>% filter(!is.na(value)) %>%
  mutate(months = log(as.numeric(substr(variable,6,8))),  #log of months
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve = log(1 - (value/100)), #log of 1-VE
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "2-dose", "3-dose"), ordered = TRUE),
         age_group = Age, #factor(Age, levels = c("18-49 years", "50-64 years", "65+ years")), #CHANGE!!!!
         study = as.factor(Study)) 

#Linear relationship check
attach(long_data)
plot(months, ve)

#create separate cases and severe outcomes dfs
ve_severe <- long_data_lower %>% filter(Outcome == "Severe Outcomes")

#Linear mixed effects model
severe_model <- lmer(ve ~ months + num_doses + age_group + prior_inf + (months|study), data = ve_severe)
summary(severe_model)

ranef(severe_model)

#Prediction for old model 
new_data <- ve_severe %>% filter(variable == "Month3") %>% dplyr::select('prior_inf', 'age_group', 'num_doses', 'study')
# adding_data <- new_data[16:18,] %>% mutate(prior_inf = 1,
#                                            study = "Unknown")
# prediction_data <- rbind(new_data, adding_data) %>% filter(num_doses %in% c("3-dose", "4+ doses"))
prediction_data <- new_data[rep(seq_len(nrow(new_data)), 25), ]

prediction_data$months <- rep(log(c(0.5, 1:24)), each = 15)

#Prediction for non-agestratigied new model
# prediction_data <- data.frame(prior_inf = rep(c(0:1), 25),
#                               study = rep(c("Bobrovitz et al (Lancet, 2023); Ferdinands et al (BMJ, 2022)", "Bobrovitz et al (Lancet, 2023)"), 25),
#                               months = rep(log(c(0.5, 1:24)), each = 2))

#Prediction for age-strat new model
# prediction_data <- data.frame(prior_inf = rep(c(0,0,0,1), 25),
#                               study = rep(c("Bobrovitz et al (Lancet, 2023); Ferdinands et al (BMJ, 2022)","Bobrovitz et al (Lancet, 2023); Ferdinands et al (BMJ, 2022)","Bobrovitz et al (Lancet, 2023); Ferdinands et al (BMJ, 2022)", "Bobrovitz et al (Lancet, 2023)"), 25),
#                               months = rep(log(c(0.5, 1:24)), each = 4),
#                               age_group = rep(c("18-49 years", "50-64 years", "65+ years", "all ages"), 25))

preds <- predict(severe_model, newdata = prediction_data, allow.new.levels = TRUE)

prediction_data$ve_pred <- 1 - exp(preds)

#Plotting predicted waning data
plot_data <- prediction_data %>% 
  mutate(group = ifelse(prior_inf == 1, 
                        paste0(age_group,", 3-dose, Prior infection, mean"),#"3-doses; Prior infection", #paste0(num_doses, ", ", "Prior infection, "),
                        paste0(age_group, ", 3-dose No prior infection, mean")),#paste0("3-doses; No prior infection\n", age_group)), 
         ve_pred = ve_pred * 100)  %>%
  filter(num_doses == '3-dose')

#Combining with observed severe outcome ve
observed_ve <- ve_severe %>% dplyr::select(prior_inf,age_group, num_doses, study, months, value) %>%
  #filter(num_doses == '3-dose') %>%
  rename(ve_pred = value) %>%
  mutate(group = ifelse(prior_inf == 1, 
                        paste0('Observed Data; \n', num_doses, ", ", age_group, ", Prior infection"), 
                        paste0("Observed Data; \n", num_doses, ", ", age_group, ", No prior infection\n")))

combined_ve <- rbind(plot_data, observed_ve)


combined_95UI <- rbind(plot_data_lower, plot_data_upper)
combined_all <- rbind(combined_95UI, plot_data)
#Combining with bivalent relative coverage
absolute_ve <- function(relative_ve, absolute_ve) {
  relative_rr <- 1 - relative_ve
  absolute_rr_lower <- 1 - absolute_ve
  absolute_rr_upper <- relative_rr * absolute_rr_lower
  absolute_ve_upper <- 1 - absolute_rr_upper
  return(absolute_ve_upper)
}

bivalent_data <- merge(plot_data, waning_relative_4dose, by = "months", all.x = TRUE) %>% mutate(ve_biv = absolute_ve(relative_ve/100, ve_pred/100)) %>% 
  dplyr::select(prior_inf, age_group, num_doses, study, months, ve_biv) %>% rename(ve_pred = ve_biv) %>%
  mutate(group = ifelse(prior_inf == 1, 
                        paste0(age_group, ", 4-dose; Prior infection"), 
                        paste0(age_group, ", 4-dose; No prior infection\n")),
         ve_pred = ve_pred * 100)

combined_ve <- rbind(plot_data, bivalent_data)

combined_95UI <- rbind(plot_data_upper, plot_data_lower)

ggplot(combined_all %>% filter(age_group == '65+ years'), aes(x = exp(months), y = ve_pred, color = group)) +
  geom_line(size = 1) +
  ggtitle(paste0("Protective Effectiveness Waning (95% UI) \nOutcomes: Severe Outcomes, Age: 65+ years\n Fitting model on 0-3 doses")) +
  ylab("Protective Effectiveness (%)")+
  ylim(0, 100) +
  scale_x_continuous(name="Months since vaccination/infections", limits=c(0, 25), breaks = c(0:24)) +
  # scale_color_manual(values = c(
  #   "dodgerblue2", "#E31A1C", # red
  #   "green4",
  #   "#6A3D9A", # purple
  #   "#FF7F00", # orange
  #   "black", "gold1",
  #   "skyblue2", "#FB9A99", # lt pink
  #   "palegreen2",
  #   "#CAB2D6", # lt purple
  #   "#FDBF6F", # lt orange
  #   "gray70", "khaki2",
#   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#   "darkturquoise", "green1", "yellow4", "yellow3",
#   "darkorange4", "brown"
# ))
scale_color_manual(values = c("lightskyblue1",
                              "dodgerblue2",
                              "royalblue1",
                              "tomato",
                              "red",
                              "red4")) 


write.csv(prediction_data, "data/ve_waning_prediction_severe_95UI_lower.csv")
