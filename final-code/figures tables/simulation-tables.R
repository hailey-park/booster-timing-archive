###################################################################################################
#Title: Plots of Simulation Results
#Author: Hailey Park
#Date: April 17, 2023
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


###################################################################################################
#Plotting Severe Outcome Incidence over Time
age_18_49_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-18-49 years-mean.csv")[,-1]
age_18_49_lb_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-18-49 years-lower.csv")[,-1]
age_18_49_ub_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-18-49 years-upper.csv")[,-1]

age_50_64_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-50-64 years-mean.csv")[,-1]
age_50_64_lb_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-50-64 years-lower.csv")[,-1]
age_50_64_ub_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-50-64 years-upper.csv")[,-1]

age_65_74_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-65-74 years-mean.csv")[,-1]
age_65_74_lb_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-65-74 years-lower.csv")[,-1]
age_65_74_ub_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-65-74 years-upper.csv")[,-1]

age_75plus_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-75+ years-mean.csv")[,-1]
age_75plus_lb_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-75+ years-lower.csv")[,-1]
age_75plus_ub_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-75+ years-upper.csv")[,-1]

immunocomp_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-immuno-mean.csv")[,-1]
immunocomp_lb_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-immuno-lower.csv")[,-1]
immunocomp_ub_1booster <- read.csv("simulation-results/sensitivity/no-waning/1Booster-immuno-upper.csv")[,-1]

age_18_49_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-18-49 years-mean.csv")[,-1]
age_18_49_lb_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-18-49 years-lower.csv")[,-1]
age_18_49_ub_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-18-49 years-upper.csv")[,-1]

age_50_64_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-50-64 years-mean.csv")[,-1]
age_50_64_lb_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-50-64 years-lower.csv")[,-1]
age_50_64_ub_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-50-64 years-upper.csv")[,-1]

age_65_74_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-65-74 years-mean.csv")[,-1]
age_65_74_lb_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-65-74 years-lower.csv")[,-1]
age_65_74_ub_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-65-74 years-upper.csv")[,-1]

age_75plus_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-75+ years-mean.csv")[,-1]
age_75plus_lb_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-75+ years-lower.csv")[,-1]
age_75plus_ub_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-75+ years-upper.csv")[,-1]

immunocomp_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-immuno-mean.csv")[,-1]
immunocomp_lb_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-immuno-lower.csv")[,-1]
immunocomp_ub_annualbooster <- read.csv("simulation-results/sensitivity/no-waning/annualBooster-immuno-upper.csv")[,-1]


age_18_49_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-18-49 years-mean.csv")[,-1]
age_18_49_lb_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-18-49 years-lower.csv")[,-1]
age_18_49_ub_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-18-49 years-upper.csv")[,-1]

age_50_64_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-50-64 years-mean.csv")[,-1]
age_50_64_lb_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-50-64 years-lower.csv")[,-1]
age_50_64_ub_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-50-64 years-upper.csv")[,-1]

age_65_74_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-65-74 years-mean.csv")[,-1]
age_65_74_lb_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-65-74 years-lower.csv")[,-1]
age_65_74_ub_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-65-74 years-upper.csv")[,-1]

age_75plus_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-75+ years-mean.csv")[,-1]
age_75plus_lb_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-75+ years-lower.csv")[,-1]
age_75plus_ub_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-75+ years-upper.csv")[,-1]

immunocomp_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-immuno-mean.csv")[,-1]
immunocomp_lb_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-immuno-lower.csv")[,-1]
immunocomp_ub_biannualbooster <- read.csv("simulation-results/sensitivity/no-waning/biannualBooster-immuno-upper.csv")[,-1]

sum(age_18_49_lb_1booster)/(2*1000000) *100000
sum(age_18_49_ub_1booster)/(2*1000000)*100000
sum(age_50_64_lb_1booster)/(2*1000000) *100000
sum(age_50_64_ub_1booster)/(2*1000000)*100000
sum(age_65_74_lb_1booster)/(2*1000000) *100000
sum(age_65_74_ub_1booster)/(2*1000000)*100000
sum(age_75plus_lb_1booster)/(2*1000000) *100000
sum(age_75plus_ub_1booster)/(2*1000000)*100000
sum(immunocomp_lb_1booster)/(2*1000000) *100000
sum(immunocomp_ub_1booster)/(2*1000000)*100000

sum(age_18_49_lb_annualbooster)/(2*1000000) *100000
sum(age_18_49_ub_annualbooster)/(2*1000000)*100000
sum(age_50_64_lb_annualbooster)/(2*1000000) *100000
sum(age_50_64_ub_annualbooster)/(2*1000000)*100000
sum(age_65_74_lb_annualbooster)/(2*1000000) *100000
sum(age_65_74_ub_annualbooster)/(2*1000000)*100000
sum(age_75plus_lb_annualbooster)/(2*1000000) *100000
sum(age_75plus_ub_annualbooster)/(2*1000000)*100000
sum(immunocomp_lb_annualbooster)/(2*1000000) *100000
sum(immunocomp_ub_annualbooster)/(2*1000000)*100000

sum(age_18_49_lb_biannualbooster)/(2*1000000) *100000
sum(age_18_49_ub_biannualbooster)/(2*1000000)*100000
sum(age_50_64_lb_biannualbooster)/(2*1000000) *100000
sum(age_50_64_ub_biannualbooster)/(2*1000000)*100000
sum(age_65_74_lb_biannualbooster)/(2*1000000) *100000
sum(age_65_74_ub_biannualbooster)/(2*1000000)*100000
sum(age_75plus_lb_biannualbooster)/(2*1000000) *100000
sum(age_75plus_ub_biannualbooster)/(2*1000000)*100000
sum(immunocomp_lb_biannualbooster)/(2*1000000) *100000
sum(immunocomp_ub_biannualbooster)/(2*1000000)*100000

#Tables
#NewResults: Bivalent booster effect excluded
# inspection <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years", "Immunocompromised"),
#                          total_1booster = c(2158,4365,11491,30861,24194),
#                          total_annual_booster = c(1751,3518,9093,24638,21839),
#                          total_biannual_booster = c(1417,2860,7041,19296,19785))
# 


#Tables
age_18_49_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(sum(age_18_49_1booster),
                                                     sum(age_18_49_annualbooster),
                                                     sum(age_18_49_biannualbooster))) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))


age_50_64_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(sum(age_50_64_1booster),
                                                     sum(age_50_64_annualbooster),
                                                     sum(age_50_64_biannualbooster))) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

age_65_74_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(sum(age_65_74_1booster),
                                                     sum(age_65_74_annualbooster),
                                                     sum(age_65_74_biannualbooster))) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

age_75plus_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                               total_severe_covid = c(sum(age_75plus_1booster),
                                                      sum(age_75plus_annualbooster),
                                                      sum(age_75plus_biannualbooster))) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

immunocompromised_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
  total_severe_covid = c(sum(immunocomp_1booster),
                         sum(immunocomp_annualbooster),
                         sum(immunocomp_biannualbooster))) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))


