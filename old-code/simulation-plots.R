###################################################################################################
#Title: Plots of Simulation Results
#Author: Hailey Park
#Date: April 17, 2023
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


###################################################################################################
#Plotting Severe Outcome Incidence over Time
age_18_49 <- read.csv("simulation-results/final/main/noBooster-18-49 years-mean.csv")[,-1]
age_18_49_lb <- read.csv("simulation-results/final/main/noBooster-18-49 years-lower.csv")[,-1]
age_18_49_ub <- read.csv("simulation-results/final/main/noBooster-18-49 years-upper.csv")[,-1]

age_50_64 <- read.csv("simulation-results/final/main/noBooster-50-64 years-mean.csv")[,-1]
age_50_64_lb <- read.csv("simulation-results/final/main/noBooster-50-64 years-lower.csv")[,-1]
age_50_64_ub <- read.csv("simulation-results/final/main/noBooster-50-64 years-upper.csv")[,-1]

age_65_74 <- read.csv("simulation-results/final/main/noBooster-65-74 years-mean.csv")[,-1]
age_65_74_lb <- read.csv("simulation-results/final/main/noBooster-65-74 years-lower.csv")[,-1]
age_65_74_ub <- read.csv("simulation-results/final/main/noBooster-65-74 years-upper.csv")[,-1]

age_75plus <- read.csv("simulation-results/final/main/noBooster-75+ years-mean.csv")[,-1]
age_75plus_lb <- read.csv("simulation-results/final/main/noBooster-75+ years-lower.csv")[,-1]
age_75plus_ub <- read.csv("simulation-results/final/main/noBooster-75+ years-upper.csv")[,-1]

immunocomp <- read.csv("simulation-results/final/main/noBooster-immuno-mean.csv")[,-1]
immunocomp_lb <- read.csv("simulation-results/final/main/noBooster-immuno-lower.csv")[,-1]
immunocomp_ub <- read.csv("simulation-results/final/main/noBooster-immuno-upper.csv")[,-1]

sum(colSums((age_18_49 %>% filter(prior_inf == 1))[, (7:31)]))
sum(colSums((age_18_49 %>% filter(prior_inf == 0))[, (7:31)]))

sum(colSums((age_50_64 %>% filter(prior_inf == 1))[, (7:31)]))
sum(colSums((age_50_64 %>% filter(prior_inf == 0))[, (7:31)]))

sum(colSums((age_65_74 %>% filter(prior_inf == 1))[, (7:31)]))
sum(colSums((age_65_74 %>% filter(prior_inf == 0))[, (7:31)]))

sum(colSums((age_75plus %>% filter(prior_inf == 1))[, (7:31)]))
sum(colSums((age_75plus %>% filter(prior_inf == 0))[, (7:31)]))

sum(colSums((immunocomp %>% filter(prior_inf == 1))[, (7:31)]))
sum(colSums((immunocomp %>% filter(prior_inf == 0))[, (7:31)]))


age_18_49_95ui <- data.frame(
                              mean = colSums(age_18_49[, (7:31)]),
                              lower = (colSums(age_18_49_lb[, (7:31)])),
                              upper = (colSums(age_18_49_ub[, (7:31)])),
                             # mean = age_18_49,
                             # lower = age_18_49_lb,
                             # upper = age_18_49_ub,
                             age_group = "18-49 years")

age_50_64_95ui <- data.frame(
                             mean = colSums(age_50_64[, (7:31)]),
                             lower = colSums(age_50_64_lb[, (7:31)]),
                             upper = colSums(age_50_64_ub[, (7:31)]),
                             # mean = age_50_64,
                             # lower = age_50_64_lb,
                             # upper = age_50_64_ub,
                             age_group = "50-64 years")

age_65_74_95ui <- data.frame(
                               mean = colSums(age_65_74[, (7:31)]),
                               lower = colSums(age_65_74_lb[, (7:31)]),
                               upper = colSums(age_65_74_ub[, (7:31)]),
  # mean = age_65_74,
  # lower = age_65_74_lb,
  # upper = age_65_74_ub,
                             age_group = "65-74 years")

age_75plus_95ui <- data.frame(
                                mean = colSums(age_75plus[, (7:31)]),
                                lower = colSums(age_75plus_lb[, (7:31)]),
                                upper = colSums(age_75plus_ub[, (7:31)]),
                                # mean = age_75plus,
                                # lower = age_75plus_lb,
                                # upper = age_75plus_ub,
                              age_group = "75+ years")

immunocomp_95ui <- data.frame(
                              mean = colSums(immunocomp[, (7:31)]),
                              lower = colSums(immunocomp_lb[, (7:31)]),
                              upper = colSums(immunocomp_ub[, (7:31)]),
                               # mean = immunocomp,
                               # lower = immunocomp_lb,
                               # upper = immunocomp_ub,
                             age_group = "Immunocompromised")


combined <- rbind(age_18_49_95ui, age_50_64_95ui, age_65_74_95ui, age_75plus_95ui, immunocomp_95ui)
combined$months <- rep((0:24), 5)

plot_data <- combined %>% mutate(mean = mean/10, lower = lower/10, upper = upper/10)

ggplot(plot_data, aes(months)) + 
  geom_line(aes(y=mean, color = age_group), size = 1.5) + 
  #geom_line(aes(y = upper, color = age_group), linetype = "dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = age_group), alpha=0.4) +
  xlab("Time (months)") + 
  ylab("Non-Severe COVID-19 Incidence") +
  labs(color='Risk Group') +
  guides(fill = "none") +
  ylim(0,1700)+
  #ylim(0,190) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=18))
  #ggtitle("Seroprevalence: 0% \nNo Non-Severe Infections")


###################################################################################################
#Tables
#NewResults: Bivalent booster effect excluded
inspection <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years", "Immunocompromised"),
                         total_1booster = c(1963,4348,11093,29674,22569),
                         total_annual_booster = c(1533,3386,8697,23665,20831),
                         total_biannual_booster = c(1298,2725,6819,18621,18577))

# inspection_lower <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
#                                total_1booster = c(3212,6427,15310,41641),
#                                total_annual_booster = c(2710,5369,12523,33744),
#                                total_biannual_booster = c(2271,4292,9596,26580))
# 
# inspection_upper <- data.frame(age_group = c("18-49 years", "50-64 years", "65-74 years", "75+ years"),
#                                total_1booster = c(1253,3253,9004,24272),
#                                total_annual_booster = c(1008,2478,6637,18625),
#                                total_biannual_booster = c(765,1821,4765,12966))

sum(colSums(age_18_49[, (7:31)]))
sum(colSums(age_50_64[, (7:31)]))
sum(colSums(age_65_74[, (7:31)]))
sum(colSums(age_75plus[, (7:31)]))
sum(colSums(immunocomp[, (7:31)]))

sum(colSums(age_18_49[, (7:31)]))/(2*1000000) *100000
sum(colSums(age_50_64[, (7:31)]))/(2*1000000) *100000
sum(colSums(age_65_74[, (7:31)]))/(2*1000000) *100000
sum(colSums(age_75plus[, (7:31)]))/(2*1000000) *100000
sum(colSums(immunocomp[, (7:31)]))/(2*1000000) *100000



sum(colSums(age_18_49_lb[, (7:31)]))/(2*1000000) *100000
sum(colSums(age_50_64_lb[, (7:31)]))/(2*1000000)*100000
sum(colSums(age_65_74_lb[, (7:31)]))/(2*1000000)*100000
sum(colSums(age_75plus_lb[, (7:31)]))/(2*1000000)*100000
sum(colSums(immunocomp_lb[, (7:31)]))/(2*1000000)*100000

sum(colSums(age_18_49_ub[, (7:31)]))/(2*1000000)*100000
sum(colSums(age_50_64_ub[, (7:31)]))/(2*1000000)*100000
sum(colSums(age_65_74_ub[, (7:31)]))/(2*1000000)*100000
sum(colSums(age_75plus_ub[, (7:31)]))/(2*1000000)*100000
sum(colSums(immunocomp_ub[, (7:31)]))/(2*1000000)*100000

#Tables
age_18_49_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(1963,1533,1298)) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))


age_50_64_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(4348,3386,2725)) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

age_65_74_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                              total_severe_covid = c(11093,8697,6819)) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

age_75plus_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                               total_severe_covid = c(29674,23665,18621)) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))

immunocompromised_table <- data.frame(intervention = c("1 Booster (1 dose)", "Annual Booster (2 doses)", "Biannual Booster (4 doses)"),
                               #total_severe_covid = c(18925,17772,16138)) %>%
  total_severe_covid = c(22569,20831,18577)) %>%
  mutate(absolute_annual_risk = total_severe_covid/(2*1000000),
         ARR = max(absolute_annual_risk) - absolute_annual_risk,
         RRR = ARR/max(absolute_annual_risk),
         NNT = ceiling(rep(1000000, 3)/(max(total_severe_covid) - total_severe_covid)))


###################################################################################################
#Plotting sensitivity analyses

oneBooster <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/1Booster-immuno-mean.csv")[,-1]
oneBooster_lower <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/1Booster-immuno-lower.csv")[,-1]
oneBooster_upper <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/1Booster-immuno-upper.csv")[,-1]
annualBooster <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/annualBooster-immuno-mean.csv")[,-1]
annualBooster_lower <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/annualBooster-immuno-lower.csv")[,-1]
annualBooster_upper <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/annualBooster-immuno-upper.csv")[,-1]
biannualBooster <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/biannualBooster-immuno-mean.csv")[,-1]
biannualBooster_lower <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/biannualBooster-immuno-lower.csv")[,-1]
biannualBooster_upper <- read.csv("simulation-results/final/sensitivity/pessimistic-ve/biannualBooster-immuno-upper.csv")[,-1]


group_sensitivity<- data.frame(total_severe = c(sum(oneBooster), 
                                           sum(annualBooster),
                                           sum(biannualBooster)),
                            total_severe_lb = c(sum(oneBooster_lower), 
                                                sum(annualBooster_lower),
                                                sum(biannualBooster_lower)),
                            total_severe_ub = c(sum(oneBooster_upper), 
                                                sum(annualBooster_upper),
                                                sum(biannualBooster_upper)),
                          intervention = c("1 Booster", "Annual Booster", "Biannual Booster"),
                          sensitivity_analysis = "Pessimistic VE") %>%
  mutate(absolute_risk = total_severe/(2*1000000),
         absolute_risk_lb = total_severe_lb/(2*1000000),
         absolute_risk_ub = total_severe_ub/(2*1000000))

write.csv(group_sensitivity, "simulation-results/final/sensitivity/sensitivity-summarised/immuno-pessimistic_ve.csv")




optimistic_ve <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-optimistic_ve.csv")[,-1]
optimistic_waning <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-optimistic_waning.csv")[, -1]
pessimistic_ve <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-pessimistic_ve.csv")[,-1]
pessimistic_waning <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-pessimistic_waning.csv")[, -1]
high_inc <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-high_incidence.csv")[,-1]
low_inc <- read.csv("simulation-results/final/sensitivity/sensitivity-summarised/immuno-low_incidence.csv")[, -1]

combined <- rbind(high_inc, low_inc, optimistic_ve, optimistic_waning, pessimistic_ve, pessimistic_waning) %>% select(-c("total_severe", "total_severe_lb", "total_severe_ub"))

ggplot(combined,aes(x = intervention, y = absolute_risk * 100, color = intervention)) + 
  facet_grid(. ~ sensitivity_analysis, switch = "both") + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=absolute_risk_lb * 100, ymax=absolute_risk_ub * 100), width=.1) +
  ylab("Annual Risk of Severe COVID-19 (%)") +
  xlab("Sensitivity Analysis")+
  labs(color = "Intervention") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=12),
        axis.text.x=element_blank()) +
  ylim(0, 3)
  #ggtitle("Age Group: 18-49 years")
###################################################################################################
#Plotting Time Since


age_18_49 <- read.csv("calibration-results/adj-calibration-1mil-18-49 years.csv")[,-1]
age_50_64 <- read.csv("calibration-results/adj-calibration-1mil-50-64 years.csv")[,-1]
age_65_74 <- read.csv("calibration-results/adj-calibration-1mil-65-74 years.csv")[,-1]
age_75plus <- read.csv("calibration-results/adj-calibration-1mil-75+ years.csv")[,-1]


time_since_18_49 <- age_18_49 %>% group_by(time_since_last_dose) %>% summarise(total = n()) %>% 
  mutate(time_since_last = as.Date(time_since_last_dose),
         age_group = "18-49 years")

time_since_50_64 <- age_50_64 %>% group_by(time_since_last_dose) %>% summarise(total = n()) %>% 
  mutate(time_since_last = as.Date(time_since_last_dose),
         age_group = "50-64 years;\n65-74 years;\n75+ years")

time_since_65_74 <- age_65_74 %>% group_by(time_since_last_dose) %>% summarise(total = n()) %>% 
  mutate(time_since_last = as.Date(time_since_last_dose),
         age_group = "65-74 years; \n75+ years")

time_since_75plus <- age_75plus %>% group_by(time_since_last_dose) %>% summarise(total = n()) %>% 
  mutate(time_since_last = as.Date(time_since_last_dose),
         age_group = "75+ years")

combined <- rbind(time_since_18_49, time_since_50_64) #, time_since_65_74) #, time_since_75plus)


#plot simulated prior dose
ggplot(data = combined %>% filter(!is.na(time_since_last)), aes(x = as.Date(time_since_last), y = total, color = age_group)) +
  geom_line(size = 1.5) +
  ylab("Total (individuals)") +
  xlab("Time (months)") + 
  labs(color='Age Group') +
  scale_x_date(labels = date_format("%b %Y"),
               limits = as.Date(c('2020-09-01', '2022-08-01'))) +
  scale_y_continuous(labels = label_comma()) +
  guides(fill = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=15)) 

###################################################################################################

#Comparing 2 versions of optimistic/pessimistic waning

optimistic_1Booster <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-18-49 years-pessimistic.csv")[,-1]
###################################################################################################optimistic_annualBooster <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-18-49 years-pessimistic_v2.csv")[,-1]

pessimistic_1Booster <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-50-64 years-pessimistic.csv")[,-1]
pessimistic_annualBooster <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-50-64 years-pessimistic_v2.csv")[,-1]


optimistic_1Booster_v2 <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-65-74 years-pessimistic.csv")[,-1]
optimistic_annualBooster_v2 <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-65-74 years-pessimistic_v2.csv")[,-1]

pessimistic_1Booster_v2 <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-75+ years-pessimistic.csv")[,-1]
pessimistic_annualBooster_v2 <- read.csv("simulation-results/third-pass-1mil-sensitivity/1Booster-75+ years-pessimistic_v2.csv")[,-1]


age_18 <- data.frame(oneBooster = colSums(optimistic_1Booster[, (7:31)]),
                     oneBooster_v2 = colSums(optimistic_annualBooster[, (7:31)]),
                     age = "18-49 years")

age_50 <- data.frame(oneBooster = colSums(pessimistic_1Booster[, (7:31)]),
                     oneBooster_v2 = colSums(pessimistic_annualBooster[, (7:31)]),
                     age = "50-64 years")

age_65 <- data.frame(oneBooster = colSums(optimistic_1Booster_v2[, (7:31)]),
                     oneBooster_v2 = colSums(optimistic_annualBooster_v2[, (7:31)]),
                     age = "65-74 years")

age_75 <- data.frame(oneBooster = colSums(pessimistic_1Booster_v2[, (7:31)]),
                     oneBooster_v2 = colSums(pessimistic_annualBooster_v2[, (7:31)]),
                     age = "75+ years")

combined <- rbind(age_18, age_50, age_65, age_75) 
combined$months <- rep((0:24), 4)

plot_data <- combined %>% mutate(oneBooster = oneBooster/10, oneBooster_v2 = oneBooster_v2/10)

ggplot(plot_data, aes(months)) + 
  geom_line(aes(y=oneBooster, color = age), size = 1.5) + 
  xlab("Time (months)") + 
  ylab("Severe COVID-19 Incidence") +
  labs(color='Age Group') +
  guides(fill = "none") +
  ylim(0,250)+
  #ylim(0,490) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=18)) +
ggtitle("Pessimistic Waning (Version 1) \n Intervention: 1 Booster")
###################################################################################################
waning_data_clean <- read.csv("data/ve_waning_predictions_mean_immunocompromised.csv")[,-1] %>%
  mutate(ve_pred = ve_pred + 0.14,
         group = if_else(prior_inf == 1, "Original: Prior Infection",
                         "Original: No Prior Infection"))

waning_data_clean_pess <- waning_data_clean %>%
  mutate(#ve_pred = ve_pred - (months * 0.00416), 
    ve_pred = ve_pred - 0.1, 
         group = if_else(prior_inf == 1, "Pessimistic: Prior Infection",
                                                         "Pessimistic: No Prior Infection"))

waning_data_clean_opt <- waning_data_clean %>%
  mutate(#ve_pred = if_else(ve_pred + (months * 0.00416) > max(ve_pred), max(ve_pred), ve_pred + (months * 0.00416)), 
    ve_pred = if_else(ve_pred + 0.1 > 1, 1, ve_pred + 0.1),
         group = if_else(prior_inf == 1, "Optimistic: Prior Infection",
                         "Optimistic: No Prior Infection"))

combined <- rbind(waning_data_clean, waning_data_clean_pess, waning_data_clean_opt)

ggplot(data = combined, aes(x = months, y = ve_pred * 100, color = as.factor(group))) +
  geom_line(size = 1.5) +
  xlab("Time (months)") +
  ylab("Protective Effectiveness (%)") +
  labs(color='Waning Curve') +
  ylim(0, 100) +
  #ggtitle("Version 2: Increase/Decrease Rate of Waning")+
  scale_color_manual(values = c("goldenrod1",
                                'goldenrod1',
                                "tomato",
                                "tomato",
                                "dodgerblue1",
                                "dodgerblue1")) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14)) 




###################################################################################################
hello <- age_75plus %>% rowwise() %>% mutate(total_nonsevere = sum(c_across(32:56), na.rm = T)) 
multiple_inf <- hello %>% filter(total_nonsevere > 0 & (total_hosps + total_deaths) > 0)

###################################################################################################

immuno_waning_check <- read.csv("data/immunocompromised waning.csv")

plot_data <- melt(immuno_waning_check, id = "immunocompromised") %>%
  mutate(months = rep(0:12, each = 2))

ggplot(data = plot_data, aes(x = months, y = value, color = as.factor(immunocompromised))) +
  geom_line(size = 1.5) +
  xlab("Time (months)") +
  ylab("Protective Effectiveness (%)") +
  labs(color='Immunocompromised') +
  ylim(0, 100) +
  ggtitle("Comparison of 3-dose Waning between Immunocompetent and Immunocompromised\n Literature: Ferdinands et al (BMJ, 2022)")+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=14)) 
  


