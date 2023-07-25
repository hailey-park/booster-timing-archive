###################################################################################################
#Title: Calibration
#Author: Hailey Park
#Date: April 18, 2023
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


#Load data
cases_by_month <- read.csv("data/cases_by_month.csv")[,-1]
four_doses_by_month <- read.csv("data/four_doses_by_month.csv")[,-1]
three_doses_by_month <- read.csv("data/three_doses_by_month.csv")[,-1]
waning_data_clean <- read.csv("data/ve_waning_predictions_mean.csv")[,-1]
waning_data_clean_95UI <- read.csv("data/ve_waning_predictions_95UI.csv")[,-1]

# waning_65 <- waning_data_clean %>% filter(age_group == "65+ years", ui == "mean") %>% mutate(age_group = "65-74 years")
# waning_75 <- waning_data_clean %>% filter(age_group == "65+ years", ui == "mean") %>% mutate(age_group = "75+ years")
# 
# new_waning_data_clean <- rbind(waning_data_clean %>% filter(age_group != "65+ years", ui == "mean"), waning_65, waning_75) %>% arrange(months)
# new_waning_data_clean_lower <- rbind(waning_data_clean_95UI %>% filter(age_group != "65+ years", ui == "lower"), waning_65, waning_75) %>% arrange(months)
# 
#  write.csv(rbind(new_waning_data_clean_lower, new_waning_data_clean_upper), "data/ve_waning_predictions_95UI.csv")
#  write.csv(new_waning_data_clean, "data/ve_waning_predictions_mean.csv")


#waning_data_clean <- waning_data_clean_95UI %>% filter(ui == "lower") %>% dplyr::select(-c("ui"))

avg_incidence_adj <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                avg_inc = c(.00008, .00016, .00041, .00113))


#Clean data
four_doses_by_month$month <- as.character(four_doses_by_month$month)
three_doses_by_month$month <- as.character(three_doses_by_month$month)

#Create matrices for each age group
age_18_49_cal <- data.frame(individual = c(1:1000000),
                        age_group = '18-49 years',
                        num_doses = '3-dose',
                        prior_inf = rbinom(1000000, 1, 0.8238))

age_50_64_cal <- data.frame(individual = c(1:1000000),
                        age_group = '50-64 years',
                        num_doses = sample(c('3-dose', '4-dose'), 1000000, prob = c(0.6, 0.4), replace = TRUE),
                        
                        prior_inf = rbinom(1000000, 1, 0.6579))

age_65_74_cal <- data.frame(individual = c(1:1000000),
                        age_group = '65-74 years',
                        num_doses = sample(c('3-dose', '4-dose'), 1000000, prob = c(0.6, 0.4), replace = TRUE),
                        
                        prior_inf = rbinom(1000000, 1, 0.4681))

age_75_plus_cal <- data.frame(individual = c(1:1000000),
                          age_group = '75+ years',
                          num_doses = as.character(sample(c('3-dose', '4-dose'), 1000000, prob = c(0.6, 0.4), replace = TRUE)),
                         
                          prior_inf = rbinom(1000000, 1, 0.4681))



#Calculate time since last vaccine dose or infection
add.months= function(date,n) {seq(date, by = paste (n, "months"), length = 2)[2]}

time_since_last <- function(df) {
  #Calculate time since last dose and time since last infection
  last_dose_and_inf <- df %>% mutate(time_since_last_dose = ifelse(num_doses == "3-dose", 
                                                                   sample(three_doses_by_month$month,
                                                                          size = sum(num_doses == '3-dose'),
                                                                          prob = three_doses_by_month$perc_doses, 
                                                                          replace = TRUE),
                                                                   sample(four_doses_by_month$month, 
                                                                          size = sum(num_doses == '4-dose'),
                                                                          prob = four_doses_by_month$perc_doses, 
                                                                          replace = TRUE)),
                                     time_since_last_inf = ifelse(prior_inf == 1,
                                                                  sample(as.character(cases_by_month$month), 
                                                                         size = sum(prior_inf == 1),
                                                                         prob = cases_by_month$perc_cases,
                                                                         replace = TRUE),
                                                                  NA))
  #Reinfection for 10% of infected individuals
  reinfection <- last_dose_and_inf %>% filter(prior_inf == 1) %>%  sample_frac(.1) %>% mutate(reinf_period = interval(as.Date(time_since_last_inf),
                                                                                                                      as.Date("2022-06-01"),) %/% months(1)) %>%
    rowwise() %>% mutate(time_since_last_reinf = ifelse(reinf_period > 0,
                                                        sample(as.character(cases_by_month$month[as.Date(cases_by_month$month) >= add.months(as.Date(time_since_last_inf), 3)]),
                                                               size = 1,
                                                               prob = cases_by_month$perc_cases[as.Date(cases_by_month$month) >= add.months(as.Date(time_since_last_inf), 3)],
                                                               replace = TRUE),
                                                        NA)) %>% dplyr::select(individual, time_since_last_reinf)                                                          
  
  
  #Merge reinfection data to main df
  merged <- merge(last_dose_and_inf, reinfection, by = c("individual"), all.x = TRUE) %>% mutate(time_since_last_dose_inf = pmax(as.Date(time_since_last_dose), as.Date(time_since_last_inf), as.Date(time_since_last_reinf), na.rm =  TRUE))
  
  return(merged)
}


####################################################
#PLOT TIME SINCE LAST
age_75plus <- time_since_last(age_75_plus_cal)

time_since_last <- age_75plus %>% group_by(time_since_last_dose) %>% summarise(total = n()) %>% 
  mutate(time_since_last_dose = as.Date(time_since_last_dose))


#plot simulated prior dose
ggplot(data = time_since_last %>% filter(!is.na(time_since_last_dose)), aes(x = as.Date(time_since_last_dose), y = total)) +
  geom_line(size = 1.5) +
  ylab("Total") +
  xlab("Months") +
  #ggtitle("Simulation of Last COVID-19 Immune Event")+
  scale_x_date(labels = date_format("%b %Y")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=15)) 
####################################################
#CALIBRATION
calibration <- function(df) {
  
  #convert time since last dose/inf into months before 9/1/22
  months_since <- df %>% mutate(months_since_last_dose_inf = as.factor(interval(time_since_last_dose_inf,as.Date('2022-09-01')) %/% months(1)))
  
  #merge with protective effectiveness
  combined <- merge(merge(months_since, waning_data_clean, by.x = c("months_since_last_dose_inf", "age_group", "prior_inf"), by.y = c("months", "age_group", "prior_inf"), all.x = TRUE),
                    avg_incidence_adj, by = "age_group", all.x = TRUE)
  
  #Calibrate risk factor (Lambda)
  pe_list <- combined$ve_pred
  avg_inc <- (combined$avg_inc)[1]
  
  lambda_calibration <- function(lambda) {
    lambda_cal <- mean(lambda * (1-pe_list))
    return(((lambda_cal - avg_inc)^2))
  }
  
  lambda <- nlm(lambda_calibration, p = c(0.01))
  print(lambda$estimate)
  
  return(combined %>% mutate(lambda = lambda$estimate))
}


####################################################
#Running Calibration on Population of 1,000,000 

calibration_results <- list(age_18_49_cal, age_50_64_cal, age_65_74_cal, age_75_plus_cal) %>%
  lapply(time_since_last) %>%
  lapply(calibration)


save_results <- function(df){
  age_group <- (df$age_group)[1]
  write.csv(df, paste0("calibration-results/final/pessimistic-ve/adj-calibration-1mil-",age_group,"-mean.csv"))
}

calibration_results %>% lapply(save_results)



inspection <- calibration_results[[2]]

#Saving Lambda 
calibrated_lambda_adj <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                lambda_mean = c(0.0005555985, 0.001001004, 0.002071947, 0.005717302),
                                lambda_upper = c(0.001027599, 0.001847575,0.003823486,0.01053799),
                                lambda_lower = c(0.0003202864,0.0005684045,0.001154417,0.003182364))

calibrated_lambda_adj_highInc <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                lambda = c(0.001110927
                                           [1] 0.002001552
                                           [1] 0.004143487
                                           [1] 0.01141632))

calibrated_lambda_adj_lowInc <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                            lambda = c(0.0002773267
                                                       [1] 0.0005005139
                                                       [1] 0.001036093
                                                       [1] 0.002853599))
