###################################################################################################
#Title: Simulation with Non-Severe Infection with smoothing and immunocompromised
#Author: Hailey Park
#Date: May 10, 2023
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
waning_data_clean <- read.csv("data/optimistic-ve/ve_waning_predictions_mean_immunocompromised.csv")[,-1] 
waning_data_clean_95UI <- read.csv("data/optimistic-ve/ve_waning_predictions_95UI_immunocompromised.csv")[,-1]
hosp_death_age_stratified <- data.frame(age_group = c("immunocompromised"),
                                        perc_death = c(0.074723946))
nonsevere_infection_multipliers <- data.frame(age_group = c("immunocompromised"),
                                              multiplier = c(16.857))

#Read in matrices for each age group
immunocompromised <- read.csv("calibration-results/final/optimistic-ve/adj-calibration-1mil-immunocompromised-mean.csv")[,-1]

#clean age matrices
clean_df <- immunocompromised %>%
  dplyr::select(c("individual", "prior_inf", "months_since_last_dose_inf", "num_doses","lambda")) %>%
  mutate(age_group = "immunocompromised")

###########################################################################################
#If doing 95% UI, set ui limit
#waning_data_clean <- waning_data_clean_95UI %>% filter(ui == "lower")

# ## #optimistic/pessimistic waning
# waning_data_clean <- waning_data_clean %>%
#    mutate(ve_pred = if_else(ve_pred + (months * 0.00416) > max(ve_pred), max(ve_pred), ve_pred + (months * 0.00416)))
# #mutate(ve_pred = ve_pred - (months * 0.00416))
# #mutate(ve_pred = ve_pred - 0.1)
# #mutate(ve_pred = if_else(ve_pred + 0.1 > 1, 1, ve_pred + 0.1))

#Function for outcome occurrence based on risk (Risk = Lambda* (1 - PE))
outcome_occurrence <- function(age, inf, time, lambda, perfect_immunity_counter, death_marker) {
  pe <- rep(1, length(age))
  for (i in which(perfect_immunity_counter == 0 & death_marker == 0)) {
    pe[i] <- (waning_data_clean[waning_data_clean$age_group == age[i] & waning_data_clean$prior_inf == inf[i] & waning_data_clean$months == time[i],])$ve_pred
  }
  severe_risk <- lambda * (1 - pe)
  nonsevere_multiplier <- (nonsevere_infection_multipliers %>% filter(age_group == "immunocompromised"))$multiplier
  nonsevere_risk <- severe_risk * nonsevere_multiplier
  return(list(rbinom(length(severe_risk), 1, severe_risk), rbinom(length(nonsevere_risk), 1, nonsevere_risk)))
}


set.seed(88)
###########################################################################################

noBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == 'immunocompromised'))$perc_death

  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf 
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)

  #Iterate through each time step
  for (i in (1:25)) {
    
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
  
  write.csv(colSums(input[, (7:31)]), "simulation-results/final/sensitivity/optimistic-ve/noBooster-immuno-mean.csv")
}

list(clean_df) %>% lapply(noBoosterSimulation)

oneBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == 'immunocompromised'))$perc_death
  averted['vaccine_wave'] <- sample(c(1:3), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf 
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  
  #Iterate through each time step
  for (i in (1:25)) {

    #Staggering vaccination over 3-month window
    if(i == 1){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i == 2) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i == 3) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
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
  
  write.csv(colSums(input[, (7:31)]), "simulation-results/final/sensitivity/high-incidence/1Booster-immuno-lower.csv")
}

list(clean_df) %>% lapply(oneBoosterSimulation)


annualBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == 'immunocompromised'))$perc_death
  averted['vaccine_wave'] <- sample(c(1:3), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf 
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  
  #Iterate through each time step
  for (i in (1:25)) {
    
    #Staggering vaccination over 3-month window
    if(i %in% c(1, 13)){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i %in% c(2,14)) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i %in% c(3,15)) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
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
  
  write.csv(colSums(input[, (7:31)]), paste0("simulation-results/final/sensitivity/high-incidence/annualBooster-immuno-lower.csv"))
}

list(clean_df) %>% lapply(annualBoosterSimulation)


biannualBoosterSimulation <- function(df){
  
  #Store averted outcomes in new df
  averted <- df
  averted[sprintf("month%s",(0:24))] <- NA
  averted[sprintf("nonsevere_month%s",(0:24))] <- NA
  averted['total_deaths'] <- 0
  averted['total_hosps'] <- 0
  averted['perc_death'] <- (hosp_death_age_stratified %>% filter(age_group == 'immunocompromised'))$perc_death
  averted['vaccine_wave'] <- sample(c(1:3), nrow(averted), replace = TRUE)
  
  input<- averted
  
  #Population's info (age_group, num_doses, prior_inf, etc.) at each timestep
  age <- as.character(input$age_group)
  doses <- as.character(input$num_doses)
  inf <- input$prior_inf
  time_since_last <- input$months_since_last_dose_inf 
  lambda <- input$lambda
  prob_death <- input$perc_death
  perfect_immunity_counter <- rep(0,nrow(input)) #If non-death infection occurs, counting down perfect immunity months
  death_marker <- rep(0,nrow(input)) #If death occurs
  hosp_count <- rep(0, nrow(input))
  death_count <- rep(0, nrow(input))
  months <- c(0.5, 1:24)
  vaccine_wave <- input$vaccine_wave
  
  #Iterate through each time step
  for (i in (1:25)) {
    
    #Staggering vaccination over 3-month window
    if(i %in% c(1, 7, 13, 19)){
      first_wave_index <- which(vaccine_wave == 1)
      time_since_last[first_wave_index] <- 0
    } else if(i %in% c(2, 8, 14, 20)) {
      second_wave_index <- which(vaccine_wave == 2)
      time_since_last[second_wave_index] <- 0
    } else if (i %in% c(3, 9, 15, 21)) {
      third_wave_index <- which(vaccine_wave == 3)
      time_since_last[third_wave_index] <- 0
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
  
  write.csv(colSums(input[, (7:31)]), paste0("simulation-results/final/sensitivity/high-incidence/biannualBooster-immuno-lower.csv"))
}

list(clean_df) %>% lapply(biannualBoosterSimulation)

###########################################################################################