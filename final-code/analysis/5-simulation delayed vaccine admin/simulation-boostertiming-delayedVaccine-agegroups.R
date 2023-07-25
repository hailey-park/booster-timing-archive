###################################################################################################
#Title: Simulation with Non-Severe Infection with smoothing
#Author: Hailey Park
#Date: May 2, 2023
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
hosp_death_age_stratified <- read.csv("data/hosp_death_age_stratified_counts_adj.csv")[,-1]
nonsevere_infection_multipliers <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                                              multiplier = c(100, 39.5, 10.9, 4.3))

############################################################################
#MAKE SURE LOADING IN CORRECT WANING DATA FROM CORRECT FOLDER
waning_data_mean <- read.csv("data/ve_waning_predictions_mean.csv")[,-1] 
waning_data_95UI <- read.csv("data/ve_waning_predictions_95UI.csv")[,-1]  


#MAKE SURE YOU ARE SETTING THE CORRECT WANING CURVE FOR CALIBRATION (mean, lower, upper)
waning_data_clean <- waning_data_mean
#waning_data_clean <- waning_data_95UI %>% filter(ui == "upper")
#waning_data_clean <- waning_data_95UI %>% filter(ui == "lower")


#MAKE SURE YOU ARE READING IN THE CORRECT CALIBRATION FILE
age_18_49 <- read.csv("calibration-results/main/adj-calibration-1mil-18-49 years-mean.csv")[,-1]
age_50_64 <-  read.csv("calibration-results/main/adj-calibration-1mil-50-64 years-mean.csv")[,-1]
age_65_74 <-  read.csv("calibration-results/main/adj-calibration-1mil-65-74 years-mean.csv")[,-1]
age_75_plus <-  read.csv("calibration-results/main/adj-calibration-1mil-75+ years-mean.csv")[,-1]

############################################################################
#clean age matrices
clean_age_matrix <- function(df){
  df %>% dplyr::select(c("individual", "age_group", "prior_inf", "months_since_last_dose_inf", "num_doses","lambda"))
}

clean_df <- list(age_18_49, age_50_64, age_65_74, age_75_plus) %>%
  lapply(clean_age_matrix) 


#Function for outcome occurrence based on risk (Risk = Lambda* (1 - PE))
outcome_occurrence <- function(age, inf, time, lambda, perfect_immunity_counter, death_marker) {
  
  pe <- rep(1, length(age))
  
  # for (i in which(perfect_immunity_counter == 0 & death_marker == 0)) {
  #   pe[i] <- (waning_data_clean[waning_data_clean$age_group == age[i] & waning_data_clean$prior_inf == inf[i] & waning_data_clean$months == time[i],])$ve_pred
  # }

  #Creating a df of individuals eligible for infection to merge with waning_data_clean to get protection at specific time point
  index_individuals_eligible <- which(perfect_immunity_counter == 0 & death_marker == 0)
  df_individuals_eligible <- data.frame(index_individual = index_individuals_eligible,
                                        age_group = age[index_individuals_eligible],
                                        prior_inf = inf[index_individuals_eligible],
                                        months = time[index_individuals_eligible])

  df_protection <- merge(df_individuals_eligible, waning_data_clean %>% dplyr::select(-c("ui")), by = c("age_group", "prior_inf", "months"), all.x = TRUE) %>% arrange(index_individual)
  print(df_protection %>% filter(index_individual == 546))
  pe[index_individuals_eligible] <- df_protection$ve_pred
  
  severe_risk <- lambda * (1 - pe)
  nonsevere_multiplier <- (nonsevere_infection_multipliers %>% filter(age_group == age[1]))$multiplier
  nonsevere_risk <- severe_risk * nonsevere_multiplier
  return(list(rbinom(length(severe_risk), 1, severe_risk), rbinom(length(nonsevere_risk), 1, nonsevere_risk), mean(pe)))
}

set.seed(88)
###########################################################################################

oneBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df %>% arrange(individual)
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  age_info <- averted$age_group[1]
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == age_info))$perc_death
  averted['vaccine_wave'] <- sample(c(1:6), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  index_time_since_last_3_months <- which(time_since_last == 3) #Individuals infected in 3 months preceding start of sim have perfect immunity at start
  index_time_since_last_2_months <- which(time_since_last == 2)
  index_time_since_last_1_months <- which(time_since_last == 1)
  perfect_immunity_counter[index_time_since_last_1_months] <- 3
  perfect_immunity_counter[index_time_since_last_2_months] <- 2
  perfect_immunity_counter[index_time_since_last_3_months] <- 1
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  average_pe <- rep(0, length(months))
  
  #Iterate through each time step
  for (i in (1:25)) {
    print(i)
    #Staggering vaccination over 6-month window
    if(i == 1){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i == 2) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i == 3) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
    }else if (i == 4) {
      fourth_wave_index <- which(vaccine_wave == 4)
      time_since_last[fourth_wave_index] <- 0
    }else if (i == 5) {
      fifth_wave_index <- which(vaccine_wave == 5)
      time_since_last[fifth_wave_index] <- 0
    }else if (i == 6) {
      sixth_wave_index <- which(vaccine_wave == 6)
      time_since_last[sixth_wave_index] <- 0
    }
    
    time_since_last[time_since_last >= 24] <- 24     #Assuming that >24 month waning is same as 24 month waning pe
    
    month <- months[time_since_last + 1]
    
    #Do outcomes occur?
    outcomes <- outcome_occurrence(age, inf, month, lambda, perfect_immunity_counter, death_marker)
    severe_outcomes <- outcomes[[1]]
    nonsevere_outcomes <- outcomes[[2]]
    average_pe[i] <- outcomes[[3]]
    
    #If no outcome occurs, increase time since last
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1
    
    #Decrease 1 from perfect immunity counter (if applicable)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1
    
    #If outcome occurs, 
    #change their prior infection status to 1, time since last to 0, perfect immunity counter to 3
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1
    time_since_last[index_outcome] <- 0
    perfect_immunity_counter[index_outcome] <- 3
    
    #Then check if severe outcome is hosp vs. death
    index_severe_outcome <- which(severe_outcomes == 1)
    death_ind <- rbinom(length(index_severe_outcome), 1, prob_death[1])
    
    #If death, cut simulation for individual (death marker)
    death_ind_index <- index_severe_outcome[which(death_ind == 1)]
    death_marker[death_ind_index] <- 1
    death_count[death_ind_index] <- death_count[death_ind_index] + 1
    
    hosp_ind_index <- index_severe_outcome[which(death_ind == 0)]
    hosp_count[hosp_ind_index] <- hosp_count[hosp_ind_index] + 1
    
    #If both severe outcome and nonsevere outcome occur in same individual, remove nonsevere outcome
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0
    
    #Add data to dataframe
    input[, i + 6] <- severe_outcomes
    input[, i + 31] <- nonsevere_outcomes
    
    
    
  }
  input$total_hosps <- hosp_count
  input$total_deaths <- death_count
  input[i,,drop = FALSE]
  
  write.csv(colSums(input[, (7:31)]), paste0("simulation-results/sensitivity/delayed-vaccine-administration/1Booster-", age_info, "-mean.csv")) ####CHANGE HERE
}

set.seed(88)
clean_df %>% lapply(oneBoosterSimulation)


annualBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df %>% arrange(individual)
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  age_info <- averted$age_group[1]
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == age_info))$perc_death
  averted['vaccine_wave'] <- sample(c(1:6), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  index_time_since_last_3_months <- which(time_since_last == 3) #Individuals infected in 3 months preceding start of sim have perfect immunity at start
  index_time_since_last_2_months <- which(time_since_last == 2)
  index_time_since_last_1_months <- which(time_since_last == 1)
  perfect_immunity_counter[index_time_since_last_1_months] <- 3
  perfect_immunity_counter[index_time_since_last_2_months] <- 2
  perfect_immunity_counter[index_time_since_last_3_months] <- 1
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  
  #Iterate through each time step
  for (i in (1:25)) {
    
    #Staggering vaccination over 6-month window
    if(i %in% c(1, 13)){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i %in% c(2,14)) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i %in% c(3,15)) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
    }else if (i %in% c(4,16)) {
      fourth_wave_index <- which(vaccine_wave == 4)
      time_since_last[fourth_wave_index] <- 0
    }else if (i %in% c(5,17)) {
      fifth_wave_index <- which(vaccine_wave == 5)
      time_since_last[fifth_wave_index] <- 0
    }else if (i %in% c(6,18)) {
      sixth_wave_index <- which(vaccine_wave == 6)
      time_since_last[sixth_wave_index] <- 0
    }
    
    time_since_last[time_since_last >= 24] <- 24     #Assuming that >24 month waning is same as 24 month waning pe
    
    month <- months[time_since_last + 1]
    
    #Do outcomes occur?
    outcomes <- outcome_occurrence(age, inf, month, lambda, perfect_immunity_counter, death_marker)
    severe_outcomes <- outcomes[[1]]
    nonsevere_outcomes <- outcomes[[2]]
    
    #If no outcome occurs, increase time since last
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1
    
    #Decrease 1 from perfect immunity counter (if applicable)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1
    
    #If outcome occurs, 
    #change their prior infection status to 1, time since last to 0, perfect immunity counter to 3
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1
    time_since_last[index_outcome] <- 0
    perfect_immunity_counter[index_outcome] <- 3
    
    #Then check if severe outcome is hosp vs. death
    index_severe_outcome <- which(severe_outcomes == 1)
    death_ind <- rbinom(length(index_severe_outcome), 1, prob_death[1])
    
    #If death, cut simulation for individual (death marker)
    death_ind_index <- index_severe_outcome[which(death_ind == 1)]
    death_marker[death_ind_index] <- 1
    death_count[death_ind_index] <- death_count[death_ind_index] + 1
    
    hosp_ind_index <- index_severe_outcome[which(death_ind == 0)]
    hosp_count[hosp_ind_index] <- hosp_count[hosp_ind_index] + 1
    
    #If both severe outcome and nonsevere outcome occur in same individual, remove nonsevere outcome
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0
    
    #Add data to dataframe
    input[, i + 6] <- severe_outcomes
    input[, i + 31] <- nonsevere_outcomes
    
    
    
  }
  input$total_hosps <- hosp_count
  input$total_deaths <- death_count
  input[i,,drop = FALSE]
  
  write.csv(colSums(input[, (7:31)]), paste0("simulation-results/sensitivity/delayed-vaccine-administration/annualBooster-", age_info, "-mean.csv")) ####CHANGE HERE
}

set.seed(88)
clean_df %>% lapply(annualBoosterSimulation)


biannualBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df %>% arrange(individual)
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  age_info <- averted$age_group[1]
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == age_info))$perc_death
  averted['vaccine_wave'] <- sample(c(1:6), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  index_time_since_last_3_months <- which(time_since_last == 3) #Individuals infected in 3 months preceding start of sim have perfect immunity at start
  index_time_since_last_2_months <- which(time_since_last == 2)
  index_time_since_last_1_months <- which(time_since_last == 1)
  perfect_immunity_counter[index_time_since_last_1_months] <- 3
  perfect_immunity_counter[index_time_since_last_2_months] <- 2
  perfect_immunity_counter[index_time_since_last_3_months] <- 1
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  
  #Iterate through each time step
  for (i in (1:25)) {
    
    #Staggering vaccination over 6-month window
    if(i %in% c(1,7,13,19)){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i %in% c(2, 8, 14, 20)) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i %in% c(3, 9, 15, 21)) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
    }else if (i %in% c(4, 10, 16, 22)) {
      fourth_wave_index <- which(vaccine_wave == 4)
      time_since_last[fourth_wave_index] <- 0
    }else if (i %in% c(5, 11, 17, 23)) {
      fifth_wave_index <- which(vaccine_wave == 5)
      time_since_last[fifth_wave_index] <- 0
    }else if (i %in% c(6, 12, 18, 24)) {
      sixth_wave_index <- which(vaccine_wave == 6)
      time_since_last[sixth_wave_index] <- 0
    }
    
    time_since_last[time_since_last >= 24] <- 24     #Assuming that >24 month waning is same as 24 month waning pe
    
    month <- months[time_since_last + 1]
    
    
    #Do outcomes occur?
    outcomes <- outcome_occurrence(age, inf, month, lambda, perfect_immunity_counter, death_marker)
    severe_outcomes <- outcomes[[1]]
    nonsevere_outcomes <- outcomes[[2]]
    
    #If no outcome occurs, increase time since last
    index_no_outcome <- which(severe_outcomes == 0 & nonsevere_outcomes == 0)
    time_since_last[index_no_outcome] <- time_since_last[index_no_outcome] + 1
    
    #Decrease 1 from perfect immunity counter (if applicable)
    perfect_immunity_counter[perfect_immunity_counter > 0] <- perfect_immunity_counter[perfect_immunity_counter > 0] - 1
    
    #If outcome occurs, 
    #change their prior infection status to 1, time since last to 0, perfect immunity counter to 3
    index_outcome <- which(severe_outcomes == 1 | nonsevere_outcomes == 1)
    inf[index_outcome] <- 1
    time_since_last[index_outcome] <- 0
    perfect_immunity_counter[index_outcome] <- 3
    
    #Then check if severe outcome is hosp vs. death
    index_severe_outcome <- which(severe_outcomes == 1)
    death_ind <- rbinom(length(index_severe_outcome), 1, prob_death[1])
    
    #If death, cut simulation for individual (death marker)
    death_ind_index <- index_severe_outcome[which(death_ind == 1)]
    death_marker[death_ind_index] <- 1
    death_count[death_ind_index] <- death_count[death_ind_index] + 1
    
    hosp_ind_index <- index_severe_outcome[which(death_ind == 0)]
    hosp_count[hosp_ind_index] <- hosp_count[hosp_ind_index] + 1
    
    #If both severe outcome and nonsevere outcome occur in same individual, remove nonsevere outcome
    index_both_outcome <- which(severe_outcomes == 1 & nonsevere_outcomes == 1)
    nonsevere_outcomes[index_both_outcome] <- 0
    
    #Add data to dataframe
    input[, i + 6] <- severe_outcomes
    input[, i + 31] <- nonsevere_outcomes
    
    
    
  }
  input$total_hosps <- hosp_count
  input$total_deaths <- death_count
  input[i,,drop = FALSE]
  
  write.csv(colSums(input[, (7:31)]), paste0("simulation-results/sensitivity/delayed-vaccine-administration/biannualBooster-", age_info, "-mean.csv")) ####CHANGE HERE
}

set.seed(88)
clean_df %>% lapply(biannualBoosterSimulation)

###########################################################################################