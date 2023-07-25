###################################################################################################
#Title: Waning VE Model-Exponential Decay
#Author: Hailey Park
#Date: March 25, 2023
###################################################################################################

rm(list=ls())

setwd("~/UCSF Research/booster-timing/updated-code")

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
waning_data_95UI <- read.csv("data/waning_data_absolute_new_agestrat_alldose_95UI.csv")

#Separate lower and upper 95% UI estimates
lower_95UI <- waning_data_95UI %>% filter(VE.Type == "lower")
upper_95UI <- waning_data_95UI %>% filter(VE.Type == "upper")


#Reformat data to long and clean data
long_data <- melt(waning_data_age_strat_alldose) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  #months
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve =  value/100, #log of 1-VE
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "2-dose", "3-dose"), ordered = TRUE),
         age_group = Age, #factor(Age, levels = c("18-49 years", "50-64 years", "65+ years")), #CHANGE!!!!
         study = as.factor(References)) %>%
  dplyr::select(c(months, prior_inf, ve, num_doses, age_group, study)) %>%
  mutate(ve_adj = if_else(prior_inf == 1,
                         #(ve - (months * 0.0008)^(0.5)),
                         (ve - (months * 0.0009375)^(0.5)),
                          ve),
         month_input = log(months),
         ve_input = log(1 - ve_adj))

#Fit models
severe_model <- lmer(ve_input ~ month_input + num_doses + age_group + prior_inf + (month_input|study), data = long_data)
summary(severe_model)

#Prediction for old model 
new_data <- long_data %>% filter(months == 3) %>% dplyr::select('prior_inf', 'age_group', 'num_doses', 'study')
prediction_data <- new_data[rep(seq_len(nrow(new_data)), 25), ]
prediction_data$month_input <- rep(log(c(0.5, 1:24)), each = 15)
prediction_data$months <- rep((c(0.5, 1:24)), each = 15)


preds <- predict(severe_model, newdata = prediction_data, allow.new.levels = TRUE)

preds_95ui <- predictInterval(severe_model, 
                              newdata = prediction_data, 
                              n.sims = 99999, 
                              level = 0.95)

prediction_data$ve_pred_mean <- 1 - exp(preds_95ui$fit)
prediction_data$ve_pred_lower <- 1 - exp(preds_95ui$upr)
prediction_data$ve_pred_upper <- 1 - exp(preds_95ui$lwr)


combined <- melt(prediction_data %>% dplyr::select(-c("study", "month_input")) %>% filter(num_doses == "3-dose"),
                 id = c("months", "prior_inf", "age_group","num_doses")) %>%
  mutate(group = if_else(prior_inf == 1,
                         paste0(variable, ": Prior Infection"),
                         paste0(variable, ": No Prior Infection")),
         value = if_else(group == "ve_pred_lower: Prior Infection",
                         value + 0.02, value))

combined$group <- factor(combined$group, levels = c("ve_pred_lower: No Prior Infection",
                                           "ve_pred_mean: No Prior Infection",
                                           "ve_pred_upper: No Prior Infection",
                                           "ve_pred_lower: Prior Infection",
                                           "ve_pred_mean: Prior Infection",
                                           "ve_pred_upper: Prior Infection"))

#Plotting predicted waning data

ggplot(combined %>% filter(age_group == '65+ years'), aes(x = months, y = value * 100, color = group)) +
  geom_line(size = 1) +
  #ggtitle(paste0("Protective Effectiveness Waning (95% UI) \nOutcomes: Severe Outcomes, Age: 18-49 years\n Fitting model on 0-3 doses")) +
  ylab("Protective Effectiveness (%)")+
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  # geom_line(data = (long_data %>% filter(age_group == "18-49 years", num_doses == "3-dose") %>%
  #                     mutate(strata = if_else(prior_inf == 1,  "Literature Estimates: Prior Infection", "Literature Estimates: No Prior Infection"))), aes(x = months, y = ve * 100, color = strata), size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0))+
  scale_color_manual(values = c(
    # "purple",
    #                             "yellow",
                                "lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4"),
                     labels = c(
                       # "Literature Estimates: No Prior Infection", 
                       #          "Literature Estimates: Prior Infection",
                                "Lower 95% UI: No Prior Infection",
                                "Mean: No Prior Infection", 
                                "Upper 95% UI: No Prior Infection",
                                "Lower 95% UI: Prior Infection",
                                "Mean: Prior Infection",
                                "Upper 95% UI: Prior Infection")) 



#Write to csv
csv_data <- combined %>% filter(variable == "ve_pred_mean") %>% rename(ve_pred = value, ui = variable) %>%
  mutate(ui = if_else(ui == "ve_pred_mean", "mean", "upper")) %>%
  dplyr::select(-c("num_doses","group"))

write.csv(csv_data, "data/ve_waning_predictions_mean.csv")
###################################################################################################
#Check relative effectiveness

#Relative VE Function
relative_ve_fn <- function(baseline, upper) {
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  return(relative_ve)
}

adj_relative_waning <- combined %>% filter(age_group == "50-64 years", variable == "ve_pred_mean") %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = 14, default = 24)) %>% 
  dplyr::select(group, months, months_baseline, value)

waning_data <- combined %>% filter(age_group == "50-64 years", variable == "ve_pred_mean") %>%
  arrange(months) %>%
  rename(ve_pred_baseline = value) %>% dplyr::select(group, months, ve_pred_baseline)

adj_waning_clean <- merge(adj_relative_waning, waning_data, by.x = c("months_baseline", "group"),
                          by.y = c("months", "group"), all.x = TRUE) %>% 
  mutate(months = if_else(months < 1, 0.5, round(months)),
         relative_ve = relative_ve_fn(ve_pred_baseline, value),
         group = if_else(group == "ve_pred_mean: Prior Infection", "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  dplyr::select(months, relative_ve, group)


relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:5), 2),
                                   relative_ve = c(.6471, .4747, .4495, .4242, .3951, .3698, .771, .4241, .3872,.3424,.3016 ,.2529),
                                   prior_inf = rep(0:1, each = 6)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  dplyr::select(-prior_inf)


combined_rel <- rbind(adj_waning_clean, relative_ve_data_lit)

ggplot(combined_rel, aes(months, relative_ve*100, color = factor(group))) +
  geom_line() +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("Relative VE Comparison \nProposed Change: \nShifted Waning Rate of Hybrid Immunity data down 8%")



ve_waning_pred <- read.csv("data/ve_waning_predictions_mean_immunocompromised.csv")[,-1] %>%
  mutate(group = if_else(prior_inf == 1, "Mean: Prior Infection", "Mean: No Prior Infection")) %>%
  select(-c("age_group"))
ve_waning_pred_95UI <- read.csv("data/ve_waning_predictions_95UI_immunocompromised.csv")[,-1]%>%
  mutate(group = case_when(prior_inf == 1 & ui == "lower" ~ "Lower 95% UI: Prior Infection",
                           prior_inf == 0 & ui == "lower" ~ "Lower 95% UI: No Prior Infection",
                           prior_inf == 1 & ui == "upper" ~ "Upper 95% UI: Prior Infection",
                           TRUE ~ "Upper 95% UI: No Prior Infection")) %>%
  select(-c("age_group", "ui"))


combined <- rbind(ve_waning_pred, ve_waning_pred_95UI)
#Plotting predicted waning data

combined$group <- factor(combined$group, levels = c("Lower 95% UI: No Prior Infection",
                                                    "Mean: No Prior Infection",
                                                    "Upper 95% UI: No Prior Infection",
                                                    "Lower 95% UI: Prior Infection",
                                                    "Mean: Prior Infection",
                                                    "Upper 95% UI: Prior Infection"))

ggplot(combined, aes(x = months, y = ve_pred * 100, color = group)) +
  geom_line(size = 1) +
  ggtitle(paste0("Protective Effectiveness Waning (95% UI) \nOutcomes: Severe Outcomes, Immunocompromised\n Fitting model on 0-3 doses")) +
  ylab("Protective Effectiveness (%)")+
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0))+
  scale_color_manual(values = c("lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4")) 




























