###################################################################################################
#Title: Model Validation
#Author: Hailey Park
#Date: June 6, 2023
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
death_data  <- read.csv("data/Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Updated__Bivalent__Booster_Status.csv")
hosp_data <- read.csv("data/COVID-19Surveillance_All_Data.csv", skip = 2)
vax_data <- read.csv("data/COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv")

#Clean data
hosp_data_clean <- hosp_data %>% 
  filter(CATCHMENT == "Entire Network",
         NETWORK == "COVID-NET",
         MMWR.YEAR == 2022,
         AGE.CATEGORY %in% c("18-49 yr", "50-64 yr", "65-74 yr", "75-84 yr", "85+"),
         MMWR.WEEK %in% c(36:48)) %>% 
  mutate(month = case_when(MMWR.WEEK %in% c(36:39) ~ 1,
                           MMWR.WEEK %in% c(40:43) ~ 2,
                           MMWR.WEEK %in% c(44:48) ~ 3,
                           TRUE ~ 0),
         WEEKLY.RATE = as.numeric(WEEKLY.RATE)) %>%
  dplyr::select(AGE.CATEGORY, WEEKLY.RATE, month) %>% group_by(AGE.CATEGORY, month) %>% summarise(monthly_rate = sum(WEEKLY.RATE)) %>%
  group_by(AGE.CATEGORY) %>% summarise(avg_inc_monthly = mean(monthly_rate))


death_data_clean <- death_data %>%
  filter(outcome == "death",
         month %in% c("SEP 2022", "OCT 2022", "NOV 2022"),
         age_group %in% c("18-29", "30-49", "50-64", "65-79", "80+")) %>%
  mutate(new_age = if_else(age_group %in% c("18-29", "30-49"), "18-49", if_else(age_group %in% c("65-79", "80+"), "65+", age_group))) %>%
  filter(!is.na(vaccinated_with_outcome)) %>%
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

vax_data_clean <- vax_data %>%
  filter(Date == "11/30/2022",
         Location == "US")
