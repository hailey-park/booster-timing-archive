###################################################################################################
#Title: Calibration in Immunocompromised Population
#Author: Hailey Park
#Date: May 10, 2023
###################################################################################################

rm(list=ls())

setwd("~/UCSF Research/booster-timing/final-code")

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
avg_incidence_adj <- .0009494565097

############################################################################
#MAKE SURE LOADING IN CORRECT WANING DATA FROM CORRECT FOLDER
waning_data_mean <- read.csv("data/ve_waning_predictions_mean_immunocompromised.csv")[,-1] ##CHANGE HERE 
waning_data_95UI <- read.csv("data/ve_waning_predictions_95UI_immunocompromised.csv")[,-1]  #CHANGE HERE 


#MAKE SURE YOU ARE SETTING THE CORRECT WANING CURVE FOR CALIBRATION (mean, lower, upper)
waning_data_clean <- waning_data_mean
#waning_data_clean <- waning_data_95UI %>% filter(ui == "upper")
#waning_data_clean <- waning_data_95UI %>% filter(ui == "lower")

#Plot waning curves
waning_data_clean %>% 
  mutate(prior_inf=as.factor(prior_inf)) %>% 
  ggplot(aes(months, ve_pred, group=interaction(age_group, prior_inf), color=age_group, shape=prior_inf)) + 
  geom_line() + geom_point() +
  ylim(0,1)

############################################################################

#Clean data
four_doses_by_month$month <- as.character(four_doses_by_month$month)
three_doses_by_month$month <- as.character(three_doses_by_month$month)


#Create matrices for immunocompromised group
immunocompromised_cal <- data.frame(individual = c(1:1000000),
                                    num_doses = sample(c('3-dose', '4-dose'), 1000000, prob = c(0.6, 0.4), replace = TRUE),
                                    prior_inf = rbinom(1000000, 1, 0.6517))

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


#CALIBRATION
calibration <- function(df) {
  
  #convert time since last dose/inf into months before 9/1/22
  months_since <- df %>% mutate(months_since_last_dose_inf = as.factor(interval(time_since_last_dose_inf,as.Date('2022-09-01')) %/% months(1)))
  
  #merge with protective effectiveness
  combined <- merge(months_since, waning_data_clean, by.x = c("months_since_last_dose_inf", "prior_inf"), by.y = c("months", "prior_inf"), all.x = TRUE) %>%
    mutate(avg_inc = avg_incidence_adj)
  
  combined$ve_pred[combined$time_since_last_inf>="2022-06-01"] <- 1
  combined$ve_pred[combined$time_since_last_reinf>="2022-06-01"] <- 1
  
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

calibration_results <- list(immunocompromised_cal) %>%
  lapply(time_since_last) %>%
  lapply(calibration)


save_results <- function(df){
  write.csv(df, paste0("calibration-results/main/calibration-1mil-immunocompromised-mean.csv"))
}

calibration_results %>% lapply(save_results)

