###################################################################################################
#Title: COVID Hospitalization Data Cleaning
#Author: Hailey Park
#Date: April 24, 2023
###################################################################################################

rm(list=ls())

setwd("~/UCSF Research/booster-timing")


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

#Load in data
hosp_data <- read.csv("COVID-19Surveillance_All_Data.csv", skip = 2)
vax_data <- read.csv("data/Archive__COVID-19_Vaccination_and_Case_Trends_by_Age_Group__United_States.csv")
death_data  <- read.csv("Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Updated__Bivalent__Booster_Status.csv")

#Clean data
hosp_data_clean <- hosp_data %>% 
  filter(CATCHMENT == "Entire Network",
         NETWORK == "COVID-NET",
         MMWR.YEAR == 2022,
         AGE.CATEGORY %in% c("18-49 yr", "50-64 yr", "65-74 yr", "75-84 yr", "85+"),
         MMWR.WEEK %in% c(10:35)) %>% 
  mutate(month = case_when(MMWR.WEEK %in% c(10:13) ~ 1,
                           MMWR.WEEK %in% c(14:17) ~ 2,
                           MMWR.WEEK %in% c(18:22) ~ 3,
                           MMWR.WEEK %in% c(23:26) ~ 4,
                           MMWR.WEEK %in% c(27:30) ~ 5,
                           MMWR.WEEK %in% c(31:35) ~ 6,
                           TRUE ~ 0),
         WEEKLY.RATE = as.numeric(WEEKLY.RATE)) %>%
  dplyr::select(AGE.CATEGORY, WEEKLY.RATE, month) %>% group_by(AGE.CATEGORY, month) %>% summarise(monthly_rate = sum(WEEKLY.RATE)) %>%
  group_by(AGE.CATEGORY) %>% summarise(avg_inc_monthly = mean(monthly_rate))

death_data_clean <- death_data %>%
  filter(outcome == "death",
         vaccination_status == "vaccinated",
         month %in% c("MAR 2022", "APR 2022", "MAY 2022", "JUN 2022", "JUL 2022", "AUG 2022"),
         age_group %in% c("18-29", "30-49", "50-64", "65-79", "80+")) %>%
  mutate(new_age = if_else(age_group %in% c("18-29", "30-49"), "18-49", if_else(age_group %in% c("65-79", "80+"), "65+", age_group))) %>%
  group_by(month, mmwr_week, new_age) %>% summarise(total_unvax_outcomes = sum(unvaccinated_with_outcome),
                                         total_unvax_pop = sum(unvaccinated_population),
                                         total_vax_outcomes = sum(vaccinated_with_outcome),
                                         total_vax_pop = sum(vaccinated_population)) %>%
  group_by(month, new_age) %>% summarise(total_unvax_outcomes = sum(total_unvax_outcomes),
                                                    total_unvax_pop = min(total_unvax_pop),
                                                    total_vax_outcomes = sum(total_vax_outcomes),
                                                    total_vax_pop = max(total_vax_pop)) %>%
  mutate(inc_unvax = total_unvax_outcomes/total_unvax_pop * 1000000,
         inc_vax = total_vax_outcomes/total_vax_pop * 1000000) %>%
  group_by(new_age) %>% summarise(avg_inc_unvax = mean(inc_unvax),
                                  avg_inc_vax = mean(inc_vax)) %>%
  mutate(X_times_higher = avg_inc_unvax/avg_inc_vax) %>%
  rename(age_group = new_age)

inspect <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
                      avg_inc_everyone = c(15.2, 27.8, 57.7,161.1),
                      #total_unvax = c(23316268, 3177035, 1845500, 1216265),
                      #total_vax = c(116690304, 60363680, 35064504, 23109039),
                      total_unvax = c(42101690, 10814601, 1882804, 1240626),
                      total_vax = c(97830170, 52800699, 31148839, 20524750),
                      X_times_higher = c(2.32,3.4,9.52,9.37)) %>% mutate(total_pop = total_unvax + total_vax,
                                                        left_side = (((avg_inc_everyone/100000) * total_pop)/total_vax) * X_times_higher,
                                                        unvax_cases = left_side/((X_times_higher/total_vax) + (1/total_unvax)),
                                                        vax_cases = (unvax_cases/total_unvax) * (total_vax/X_times_higher),
                                                        avg_inc_vax = (vax_cases/total_vax) * 100000)
inspect2 <- inspect %>% dplyr::select(age_group, avg_inc_everyone, avg_inc_vax)
inspect2$avg_inc_guesstimate <- c(15, 35, 80, 130)

                      
perc_unvax = c(0.18, 0.059, 0.05, 0.05),
                      total_unvax = c(23316268, 3177035, 1845500, 1216265)) %>% mutate(avg_inc_unvax_total = 102.25,
                                                                          unvax_scale_factor = total_unvax/sum(total_unvax),
                                                                          avg_inc_unvax_byAge = avg_inc_unvax_total * unvax_scale_factor)








###################################################################################################
#Looking at Cumulative Infections over 2/2022 - 9/2022

infections <- vax_data %>% mutate(month = as.Date(Date.Administered, "%m/%d/%Y")) %>%
  filter(month > as.Date("2022-02-26") & month < as.Date("2022-09-01"),
         AgeGroupVacc %in% c("18 - 24 Years", "25 - 49 Years", "50 - 64 Years", "65+ Years"))

infections_v2 <- death_data %>%
  filter(outcome == "case",
         vaccination_status == "vaccinated",
         month %in% c("MAR 2022", "APR 2022", "MAY 2022", "JUN 2022", "JUL 2022", "AUG 2022"),
         age_group %in% c("18-29", "30-49", "50-64", "65-79", "80+")) %>%
  mutate(new_age = if_else(age_group %in% c("18-29", "30-49"), "18-49", if_else(age_group %in% c("65-79", "80+"), "65+", age_group))) %>%
  group_by(month, mmwr_week, new_age) %>% summarise(total_unvax_outcomes = sum(unvaccinated_with_outcome),
                                                    total_unvax_pop = sum(unvaccinated_population),
                                                    total_vax_outcomes = sum(vaccinated_with_outcome),
                                                    total_vax_pop = sum(vaccinated_population)) %>%
  group_by(month, new_age) %>% summarise(total_unvax_outcomes = sum(total_unvax_outcomes),
                                         total_unvax_pop = min(total_unvax_pop),
                                         total_vax_outcomes = sum(total_vax_outcomes),
                                         total_vax_pop = max(total_vax_pop)) %>%
  mutate(total_outcomes = total_unvax_outcomes + total_vax_outcomes,
         total_pop = total_unvax_pop + total_vax_pop) %>%
  group_by(new_age) %>% summarise(total_outcomes = sum(total_outcomes),
                                  total_pop = mean(total_pop)) %>%
  mutate(inc_perc = total_outcomes/total_pop * 100)
  
  
  


