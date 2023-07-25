###################################################################################################
#Title: "Time Since Last" Calculation
#Author: Hailey Park
#Date: March 27, 2023
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
library(MASS)

#Read in case data
cases <- read.csv("data/case_outcomes_2yr.csv")[,-1]
vaccines <- read.csv("data/covid19vaccinesadministeredbydemographics.csv")

#Clean
cases_by_month <- cases %>% group_by(month) %>% summarise(num_cases = sum(num_cases)) %>% mutate(month = as.Date(month),
                                                                                                 perc_cases = round(num_cases/sum(num_cases) * 100, 2))

doses_by_month <- vaccines %>%  mutate(month = floor_date(as.Date(administered_date), unit = "month"), administered_date = as.Date(administered_date)) %>% 
  filter((administered_date >= '2021-09-22' & administered_date <= '2022-08-30'), demographic_category == "Age Group", demographic_value %in% c("18-49", "50-64", "65+")) %>%
  mutate(total_booster_doses = total_doses - partially_vaccinated - fully_vaccinated) %>% group_by(month) %>% summarise(total_boosted = sum(booster_recip_count), total_booster_doses = sum(total_booster_doses)) %>%
  dplyr::select(month, total_booster_doses)

total_3_dose <- ceiling(sum(doses_by_month$total_boosted) * 0.6)
total_4_dose <- sum(doses_by_month$total_boosted) - total_3_dose


#plot 
ggplot(data = cases_by_month, aes(x = month, y = num_cases)) +
  geom_line() +
  ylab("Total COVID-19 Infections") +
  xlab("Months") +
  ggtitle("Monthly COVID-19 Cases over 2 Years\n(Sept 1, 2020 - Sept 1, 2022)")

ggplot(data = doses_by_month, aes(x = month, y = total_booster_doses)) +
  geom_line() +
  ylab("Total COVID-19 Booster Doses") +
  xlab("Months") +
  ggtitle("Monthly COVID-19 Booster Doses since introduction of 1st Booster Dose \n(Sept 22, 2021 - Sept 1, 2022)")


#simulating time since last infection 
time_inf <- sample(cases_by_month$month, 
                   size = 1000,
                   prob = cases_by_month$perc_cases/100, 
                   replace = TRUE)

inspection <- data.frame(month = time_inf) %>% group_by(month) %>% summarise(total = n())

#plot simulated prior inf
ggplot(data = inspection, aes(x = month, y = total)) +
  geom_line() +
  ylab("Total COVID-19 Infections") +
  xlab("Months") +
  ggtitle("Simulated COVID-19 Cases over 2 Years\n(Sept 1, 2020 - Sept 1, 2022)")


#simulating time since last booster dose
four_doses_by_month <- doses_by_month %>% filter(month >= as.Date("2022-03-01") & month <= as.Date("2022-06-01")) %>% mutate(perc_doses = round(total_booster_doses/sum(total_booster_doses), 3))
three_doses_by_month <- doses_by_month %>% filter(month <= as.Date("2022-03-01")) %>% mutate(perc_doses = round(total_booster_doses/sum(total_booster_doses), 3))


time_4_dose <- sample(four_doses_by_month$month, 
                      size = 1000,
                      prob = four_doses_by_month$perc_doses, 
                      replace = TRUE)
time_3_dose <- sample(three_doses_by_month$month, 
                      size = 1000,
                      prob = three_doses_by_month$perc_doses, 
                      replace = TRUE)

inspection <- data.frame(month = time_4_dose) %>% group_by(month) %>% summarise(total = n())

#plot simulated prior dose
ggplot(data = inspection, aes(x = month, y = total)) +
  geom_line() +
  ylab("Total COVID-19 Booster Doses") +
  xlab("Months") +
  ggtitle("Simulated COVID-19 Booster Doses (4-dose) since introduction of 1st Booster Dose \n(Sept 22, 2021 - Sept 1, 2022)")



#write as .csv
write.csv(four_doses_by_month, "data/four_doses_by_month.csv")





