###################################################################################################
#Title: Relative VE Waning Adjustment
#Author: Hailey Park
#Date: March 25, 2023
###################################################################################################

rm(list=ls())

setwd("~/UCSF Research/booster-timing/updated-code")

#Load libraries
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

#Read in data
waning_pred_severe <- read.csv("data/ve_waning_predictions.csv")[,-1]
waning_pred_lower <- read.csv("data/ve_waning_prediction_severe_95UI_lower.csv")[,-1]
waning_pred_upper <- read.csv("data/ve_waning_prediction_severe_95UI_upper.csv")[,-1]

#Clean df

clean_df_mean <- function(df){
  waning_baseline_data <- df %>%
    mutate(age_group = if_else(age_group %in% c("65-74 years"), "65+ years", age_group),
           months_baseline = lead(months, n = 56, default = 24)) %>%
    filter(age_group %in% c("18-49 years", "50-64 years", "65+ years")) %>%
    select(prior_inf, age_group, months, months_baseline, ve_pred)
  
  waning_data <- waning_baseline_data %>% select(-c("months_baseline")) %>% rename(ve_pred_baseline = ve_pred)
  
  waning_clean <- merge(waning_baseline_data, waning_data, by.x = c("months_baseline", "prior_inf", "age_group"),
                        by.y = c("months", "prior_inf", "age_group"), all.x = TRUE) %>% mutate(months = if_else(months < 1, 0.5, round(months)))
  
  return(waning_clean)
}


clean_df_95ui <- function(df){
  waning_baseline_data <- df %>% filter(num_doses == '3-dose') %>% 
    mutate(months = exp(months), 
           months_baseline = lead(months, n = 42, default = 24)) %>% 
    select(prior_inf, age_group, months, months_baseline, ve_pred)
  
  waning_data <- waning_baseline_data %>% select(-c("months_baseline")) %>% rename(ve_pred_baseline = ve_pred)
  waning_clean <- merge(waning_baseline_data, waning_data, by.x = c("months_baseline", "prior_inf", "age_group"),
                        by.y = c("months", "prior_inf", "age_group"), all.x = TRUE) %>% mutate(months = if_else(months < 1, 0.5, round(months)))
  
  waning_24mo <- waning_baseline_data %>% filter(months > 23.5) %>% select(-c("months_baseline", "months"))
  
  waning_fillna <- merge(waning_clean %>% filter(is.na(ve_pred_baseline)), waning_24mo, by.x = c("prior_inf", "age_group"),
                        by.y = c("prior_inf", "age_group"), all.x = TRUE) %>% 
    select(-c("ve_pred_baseline")) %>% rename(ve_pred_baseline = ve_pred.y, ve_pred = ve_pred.x) %>%
    select(c("months_baseline", "prior_inf", "age_group", "months", "ve_pred", "ve_pred_baseline"))
  
  waning_final <- rbind(waning_clean %>% filter(!is.na(ve_pred_baseline)), waning_fillna)
  
  return(waning_final)
}

clean_waning_data <- c((list(waning_pred_lower, waning_pred_upper) %>%
  lapply(clean_df_95ui) ), (list(waning_pred_severe) %>% lapply(clean_df_mean)))



#Calculate relative VEs from absolute VEs
relative_ve_fn <- function(baseline, upper) {
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  return(relative_ve)
}

relative_waning <- clean_waning_data[[3]] %>% mutate(relative_ve = if_else(prior_inf == 1,
                                                                 relative_ve_fn(ve_pred_baseline, ve_pred) , #+.17,
                                                                 relative_ve_fn(ve_pred_baseline, ve_pred) ), #+.05),
                                           group = if_else(prior_inf == 1, "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  dplyr::select(months, relative_ve, ve_pred, group)

relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:5), 2),
                                   relative_ve = c(.6471, .4747, .4495, .4242, .3951, .3698, .771, .4241, .3872,.3424,.3016 ,.2529),
                                   prior_inf = rep(0:1, each = 6)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  dplyr::select(-prior_inf)


combined <- rbind(relative_waning, relative_ve_data_lit)

ggplot(relative_waning, aes(months, ve_pred*100, color = factor(group))) +
  geom_line() +
  ylab("Absolute Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0))


###################################################################################################
#Adjusting Relative VEs and recalculating the absolute waning curves

rel_adjustment <- function(df) {
  rel <- df %>% mutate(relative_ve = if_else(prior_inf == 1,
                                                relative_ve_fn(ve_pred_baseline, ve_pred) +.17,
                                                relative_ve_fn(ve_pred_baseline, ve_pred) +.05))
  return(rel)
}

relative_waning_adjusted <- clean_waning_data %>%
  lapply(rel_adjustment) 

relative_waning_mean <- relative_waning_adjusted[[3]] 
relative_waning_lower <- relative_waning_adjusted[[1]]
relative_waning_upper <- relative_waning_adjusted[[2]]

absolute_ve <- function(relative_ve, absolute_ve) {
  relative_rr <- 1 - relative_ve
  absolute_rr_upper <- 1 - absolute_ve
  absolute_rr_lower <- absolute_rr_upper/relative_rr
  absolute_ve_lower <- 1 - absolute_rr_lower
  return(absolute_ve_lower)
}

#Absolute Waning: Mean
absolute_waning_adj <- relative_waning_mean %>% mutate(abs_ve = absolute_ve(relative_ve, ve_pred)) %>% arrange(months) %>%
  mutate( abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
          abs_ve_final = if_else(abs_ve_move_down == 999,
                                 ve_pred, abs_ve_move_down),
          group = if_else(prior_inf == 1, "Mean; Prior Infection", "Mean; No Prior Infection")) %>%
  select(months, age_group, abs_ve_final, prior_inf, group)

#Absolute Waning: Lower
absolute_waning_adj_lower <- relative_waning_lower %>% mutate(abs_ve = absolute_ve(relative_ve, ve_pred)) %>% arrange(months) %>%
  mutate( abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
          abs_ve_final = if_else(abs_ve_move_down == 999,
                                 ve_pred, abs_ve_move_down),
          group = if_else(prior_inf == 1, "Lower 95% UI; Prior Infection", "Lower 95% UI; No Prior Infection")) %>%
  select(months, age_group, abs_ve_final, prior_inf, group)

#Absolute Waning: Upper
absolute_waning_adj_upper <- relative_waning_upper %>% mutate(abs_ve = absolute_ve(relative_ve, ve_pred)) %>% arrange(months) %>%
  mutate( abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
          abs_ve_final = if_else(abs_ve_move_down == 999,
                                 ve_pred, abs_ve_move_down),
          group = if_else(prior_inf == 1, "Upper 95% UI; Prior Infection", "Upper 95% UI; No Prior Infection")) %>%
  select(months, age_group, abs_ve_final, prior_inf, group)



#write.csv(absolute_waning_adj_upper, "data/ve_waning_predictions_adj_upper.csv")

# old_abs_ve <- waning_clean %>% mutate(group = if_else(prior_inf == 1, "Old Predictions; Prior Infection", "Old Predictions; No Prior Infection"),
#                                       abs_ve_final = ve_pred) %>% select(months, abs_ve_final, group)

combined <- rbind(absolute_waning_adj, absolute_waning_adj_lower, absolute_waning_adj_upper) %>%
  mutate(strata = paste0(group, ", ", age_group))
combined$group <- factor(combined$group, levels = c("Lower 95% UI; No Prior Infection",
                                                    "Mean; No Prior Infection",
                                                    "Upper 95% UI; No Prior Infection",
                                                    "Lower 95% UI; Prior Infection",
                                                    "Mean; Prior Infection",
                                                    "Upper 95% UI; Prior Infection"))

ggplot(combined %>% filter(age_group == "65+ years"), aes(months, abs_ve_final*100, color = factor(group))) +
  geom_line() +
  ylab("Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("Adjusted Absolute Waning Curves \n Age Group: 65+ Years") +
  scale_color_manual(values = c("lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4")) 


###################################################################################################
#More Adjustment To Hybrid Immunity (Pessimistic)
# combined <- rbind(clean_waning_data[[1]] %>% mutate(group = if_else(prior_inf == 1, "Lower 95% UI; Prior Infection", "Lower 95% UI; No Prior Infection")),
#                   clean_waning_data[[2]] %>% mutate(group = if_else(prior_inf == 1, "Upper 95% UI; Prior Infection", "Upper 95% UI; No Prior Infection")),
#                   clean_waning_data[[3]] %>% mutate(group = if_else(prior_inf == 1, "Mean; Prior Infection", "Mean; No Prior Infection"))) %>% 
#   mutate(strata = paste0(group, ", ", age_group))
# 
# combined$group <- factor(combined$group, levels = c("Lower 95% UI; No Prior Infection",
#                                                     "Mean; No Prior Infection",
#                                                     "Upper 95% UI; No Prior Infection",
#                                                     "Lower 95% UI; Prior Infection",
#                                                     "Mean; Prior Infection",
#                                                     "Upper 95% UI; Prior Infection"))


adjusted_waning <- combined %>% mutate(ve_pred = if_else(prior_inf == 1,
                                                             abs_ve_final - (months * 0.000416)^(0.5),
                                                         abs_ve_final),
                                       strata = paste0(group, ", ", age_group),
                                       ve_pred = if_else(group == "Upper 95% UI; Prior Infection", ve_pred + 0.02, ve_pred))


ggplot(adjusted_waning %>% filter(age_group == "65+ years"), aes(months, ve_pred*100, color = factor(group))) +
  geom_line() +
  ylab("Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("New Adjusted Absolute Waning Curves \n Age Group: 65+ Years") +
  scale_color_manual(values = c("lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4")) 


#Literature Estimates of Absolute 3-dose waning
waning_data_age_strat_alldose <- long_data <- melt(read.csv("data/waning_data_absolute_new_agestrat_alldose.csv")) %>% filter(!is.na(value)) %>%
  mutate(months = as.numeric(substr(variable,6,8)),  #months
         prior_inf = ifelse(Prior.Infection == "Yes", 1, 0),
         ve = (value), #VE
         num_doses = factor(Vaccine.Status, levels = c("Unvaccinated", "2-dose", "3-dose"), ordered = TRUE),
         age_group = Age, #factor(Age, levels = c("18-49 years", "50-64 years", "65+ years")), #CHANGE!!!!
         study = as.factor(References),
         strata = if_else(prior_inf == 1,  "Literature Estimates: Prior Infection", "Literature Estimates: No Prior Infection")) %>%
  select(months, prior_inf, ve, num_doses, age_group, strata) %>%
  filter(num_doses == '3-dose')


#Fitting to exponential decay
p <- ggplot(adjusted_waning %>% filter(strata %in% c("Lower 95% UI; No Prior Infection, 18-49 years",
                                         "Mean; No Prior Infection, 18-49 years",
                                         "Upper 95% UI; No Prior Infection, 18-49 years",
                                         "Lower 95% UI; Prior Infection, 18-49 years",
                                         "Mean; Prior Infection, 18-49 years",
                                         "Upper 95% UI; Prior Infection, 18-49 years")) %>%
         mutate(strata = factor(strata, levels =  c("Lower 95% UI; No Prior Infection, 18-49 years",
                                                 "Mean; No Prior Infection, 18-49 years",
                                                 "Upper 95% UI; No Prior Infection, 18-49 years",
                                                 "Lower 95% UI; Prior Infection, 18-49 years",
                                                 "Mean; Prior Infection, 18-49 years",
                                                 "Upper 95% UI; Prior Infection, 18-49 years"))), aes(x = months, y = ve_pred* 100, group = strata, color = strata)) + 
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE) +
  geom_line(data = (waning_data_age_strat_alldose %>% filter(age_group == "18-49 years")), aes(x = months, y = ve, color = strata)) +
  ylab("Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("New Adjusted Absolute Waning Curves \n Age Group: 18-49 Years") +
  scale_color_manual(values = c("yellow",
                                'green',
                                "lightskyblue1",
                                "tomato",
                                "dodgerblue1",
                                "red",
                                "royalblue3",
                                "red4"),
                     labels = c("Literature Estimates: Prior Infection", 
                                "Literature Estimates: No Prior Infection",
                       "Lower 95% UI; No Prior Infection",
                       "Lower 95% UI; Prior Infection",
                                "Mean; No Prior Infection",
                       "Mean; Prior Infection",
                                "Upper 95% UI; No Prior Infection",
                                "Upper 95% UI; Prior Infection"))

p

p2 <- ggplot_build(p)

hello <- data.frame(months = (p2$data[[1]])$x,
                    ve_pred = (p2$data[[1]])$y / 100,
                    group = (p2$data[[1]])$colour) %>%
  mutate(months_round = round(months, 1),
         group = case_when(group == "lightskyblue1" ~ "Lower 95% UI; No Prior Infection",
                           group == "dodgerblue1" ~ "Mean; No Prior Infection",
                           group == "royalblue3" ~ "Upper 95% UI; No Prior Infection",
                           group == "tomato" ~ "Lower 95% UI; Prior Infection",
                           group == "red" ~ "Mean; Prior Infection",
                           group == "red4" ~ "Upper 95% UI; Prior Infection",)) %>%
  filter(months_round %in% c(0.5, 1.1, 2, 2.9, 4.1, 5, 5.9, 7, 7.9, 9.1, 10, 10.9, 12.1, 13, 13.9, 15.1, 16, 16.9, 18.1, 18.9, 20.1, 21, 21.9, 23.1, 24)) %>%
  mutate(months = if_else(months_round == 0.5, 0.5, round(months_round)))

#####################################
#OLD Fitting

exponential_decay_model <- function(df){
  model <- nls(ve_pred ~ SSasymp(months, yf, y0, log_alpha), data = df)
  #predictions <- coef(model)[[1]] + (coef(model)[[2]] - coef(model)[[1]])*exp(-exp(coef(model)[[3]]) * df$months)
  predictions <- predict(model, df$months)
  return(predictions)
}

fitted_predictions <- lapply(split(combined, combined$strata), exponential_decay_model)
  
age_group_18_49 <- melt(data.frame(no_prior_inf_lower = fitted_predictions[[1]],
                                   no_prior_inf_mean = fitted_predictions[[7]],
                                   no_prior_inf_upper = fitted_predictions[[13]],
                                   prior_inf_lower = fitted_predictions[[4]],
                                   prior_inf_mean = fitted_predictions[[10]],
                                   prior_inf_upper = fitted_predictions[[16]],
                                   months = c(0.5, 1:24)), id = "months")
age_group_50_64 <- melt(data.frame(no_prior_inf_lower = fitted_predictions[[2]],
                                   no_prior_inf_mean = fitted_predictions[[8]],
                                   no_prior_inf_upper = fitted_predictions[[14]],
                                   prior_inf_lower = fitted_predictions[[5]],
                                   prior_inf_mean = fitted_predictions[[11]],
                                   prior_inf_upper = fitted_predictions[[17]],
                                   months = c(0.5, 1:24)), id = "months")
age_group_65 <- melt(data.frame(no_prior_inf_lower = fitted_predictions[[3]],
                                   no_prior_inf_mean = fitted_predictions[[9]],
                                   no_prior_inf_upper = fitted_predictions[[15]],
                                   prior_inf_lower = fitted_predictions[[6]],
                                   prior_inf_mean = fitted_predictions[[12]],
                                   prior_inf_upper = fitted_predictions[[18]],
                                   months = c(0.5, 1:24)), id = "months")

ggplot(age_group_18_49, aes(months, value*100, color = factor(variable))) +
  geom_line() +
  ylab("Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("New Adjusted Absolute Waning Curves \n Age Group: 65+ Years") +
  scale_color_manual(values = c("lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4"),
                     labels = c("Lower 95% UI; No Prior Infection",
                                "Mean; No Prior Infection",
                                "Upper 95% UI; No Prior Infection",
                                "Lower 95% UI; Prior Infection",
                                "Mean; Prior Infection",
                                "Upper 95% UI; Prior Infection")) 
#####################################

ve_waning_pred <- adjusted_waning %>% filter(!group %in% c("Mean; No Prior Infection", "Mean; Prior Infection")) %>%
  select(-c(strata, abs_ve_final))

write.csv(ve_waning_pred, "data/ve_waning_predictions_adj_95UI.csv")

#####################################

#Compare relative VEs
adj_relative_waning <- hello %>% filter(group %in% c("Mean; No Prior Infection", "Mean; Prior Infection")) %>%
  arrange(months) %>%
  mutate(months_baseline = lead(months, n = 14, default = 24)) %>% 
  select(group, months, months_baseline, ve_pred)

waning_data <- hello %>% filter(group %in% c("Mean; No Prior Infection", "Mean; Prior Infection")) %>%
  arrange(months) %>%
  rename(ve_pred_baseline = ve_pred) %>% select(group, months, ve_pred_baseline)

adj_waning_clean <- merge(adj_relative_waning, waning_data, by.x = c("months_baseline", "group"),
                      by.y = c("months", "group"), all.x = TRUE) %>% 
  mutate(months = if_else(months < 1, 0.5, round(months)),
         relative_ve = relative_ve_fn(ve_pred_baseline, ve_pred),
         group = if_else(group == "Mean; Prior Infection", "Waning Predictions; Prior Infection", "Waning Predictions; No Prior Infection")) %>%
  select(months, relative_ve, group)


relative_ve_data_lit <- data.frame(months = rep(c(0.5, 1:5), 2),
                                   relative_ve = c(.6471, .4747, .4495, .4242, .3951, .3698, .771, .4241, .3872,.3424,.3016 ,.2529),
                                   prior_inf = rep(0:1, each = 6)) %>%
  mutate(group = if_else(prior_inf == 1, "Literature Estimates; Prior Infection", "Literature Estimates; No Prior Infection")) %>%
  select(-prior_inf)


combined <- rbind(adj_waning_clean, relative_ve_data_lit)

ggplot(combined, aes(months, relative_ve*100, color = factor(group))) +
  geom_line() +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("Relative VE Comparison \nProposed Change: \nShifted Relative VE Curves to better match literature estimates \nAdjusted Rate of Waning for Hybrid Immunity by 10% \nFitted Absolute Predictions (after adjustments) to Exp Decay Curve \n")
###################################################################################################
#Fitting exponential decay function to Lin estimates to use as relative VE

relative_ve_no_prior_inf <- combined %>% mutate(prior_inf = if_else(group %in% c("Waning Predictions; No Prior Infection", "Literature Estimates; No Prior Infection"),0, 1)) %>%
  filter(prior_inf == 0)
relative_ve_prior_inf <- combined %>% filter(group %in% c("Waning Predictions; Prior Infection", "Literature Estimates; Prior Infection")) 

fit <- nls(relative_ve ~ SSasymp(months, yf, y0, log_alpha), data = relative_ve_prior_inf)
fit2 <- nls(relative_ve ~ SSasymp(months, yf, y0, log_alpha), data = relative_ve_no_prior_inf)

# Make predictions
predictions_new <- melt(data.frame(months = c(0.5, 1:24)) %>%
  mutate(relative_pred_prior = coef(fit)[[1]] + (coef(fit)[[2]] - coef(fit)[[1]])*exp(-exp(coef(fit)[[3]]) * months),
         relative_pred_no_prior = coef(fit2)[[1]] + (coef(fit2)[[2]] - coef(fit2)[[1]])*exp(-exp(coef(fit2)[[3]]) * months)), id = "months")


ggplot()+
  geom_point(data = relative_ve_no_prior_inf, aes(months, relative_ve, color = group), size = 0.5) +
  geom_point(data = relative_ve_prior_inf, aes(months, relative_ve, color = group), size = 0.5) +
geom_line(data = predictions_new, aes(months, value, color = variable), size= 1) + ylim(0, 1) +
  ylab("Relative Protective Effectiveness (%)") +
  xlab("Time (months)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  scale_color_manual(values = c("lightskyblue1",
                                "orangered",
                                "springgreen",
                                "orange",
                                "royalblue3",
                                "red4"),
                     labels = c("Literature Estimates: No Prior Infection",
                                "Literature Estimates: Prior Infection",
                                "Fitted Exp Decay Model: No Prior Infection",
                                "Fitted Exp Decay Model: Prior Infection",
                                "Old Waning Predictions: No Prior Infection",
                                "Old Waning Predictions: Prior Infections")) +
  labs(color='Group') +
  ggtitle("New Relative Waning Predictions Fitted with Exponential Decay Model\n ")
  
###################################################################################################
#Reconverting fitted relative waning curves back to absolute waning curves

final_relative_waning <- predictions_new %>% mutate(prior_inf = if_else(variable == "relative_pred_prior", 1, 0)) %>% 
  select(-c("variable")) %>% rename(rel_ve_pred = value)

convert_absolute_waning <- merge(adjusted_waning, final_relative_waning, by = c("months", "prior_inf"), all.x = TRUE) %>%
  mutate(abs_ve = absolute_ve(rel_ve_pred, adj_ve_pred)) %>% arrange(months) 

#Absolute Waning: Mean
final_abs_mean <- convert_absolute_waning %>% filter(group %in% c("Mean; Prior Infection", "Mean; No Prior Infection")) %>%
  mutate(abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
         abs_ve_final = if_else(abs_ve_move_down == 999,
                                adj_ve_pred, abs_ve_move_down))

#Absolute Waning: Lower 95% UI
final_abs_lower <- convert_absolute_waning %>% filter(group %in% c("Lower 95% UI; Prior Infection", "Lower 95% UI; No Prior Infection")) %>%
  mutate(abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
         abs_ve_final = if_else(abs_ve_move_down == 999,
                                adj_ve_pred, abs_ve_move_down))

#Absolute Waning: Upper 95% UI
final_abs_upper <- convert_absolute_waning %>% filter(group %in% c("Upper 95% UI; Prior Infection", "Upper 95% UI; No Prior Infection")) %>%
  mutate(abs_ve_move_down = lag(abs_ve ,n = 42, default = 999),
         abs_ve_final = if_else(abs_ve_move_down == 999,
                                adj_ve_pred, abs_ve_move_down))



combined_final <- rbind(final_abs_mean, final_abs_lower, final_abs_upper)
combined_final$group <- factor(combined_final$group, levels = c("Lower 95% UI; No Prior Infection",
                                                    "Mean; No Prior Infection",
                                                    "Upper 95% UI; No Prior Infection",
                                                    "Lower 95% UI; Prior Infection",
                                                    "Mean; Prior Infection",
                                                    "Upper 95% UI; Prior Infection"))

ggplot(combined_final %>% filter(age_group == "65+ years"), aes(months, abs_ve_final*100, color = factor(group))) +
  geom_line() +
  ylab("Protective Effectiveness (%)") +
  xlab("Time (months)") +
  labs(color='Group') +
  ylim(0,100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12))+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  ggtitle("Adjusted Absolute Waning Curves \n Age Group: 65+ Years") +
  scale_color_manual(values = c("lightskyblue1",
                                "dodgerblue1",
                                "royalblue3",
                                "tomato",
                                "red",
                                "red4")) 






