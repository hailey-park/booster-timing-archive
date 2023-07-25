library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(lubridate)
library(reshape2)
library(stringr)


#####################################


relative_ve <- function(baseline, upper) {
  baseline_rr <- 1 - baseline
  upper_rr <- 1- upper
  relative_rr <- upper_rr / baseline_rr
  relative_ve <- 1 - relative_rr
  return(relative_ve)
}


absolute_ve <- function(relative_ve, absolute_ve) {
  relative_rr <- 1 - relative_ve
  absolute_rr_lower <- 1 - absolute_ve
  absolute_rr_upper <- relative_rr * absolute_rr_lower
  absolute_ve_upper <- 1 - absolute_rr_upper
  return(absolute_ve_upper)
}

absolute_ve(.23, .80)
relative_ve(.83, .90)



 #####################################
waning_data <- read.csv("waning_data_absolute.csv") 
  #select(-c(18:31))

plot_list = list()
counter <- 1
for (outcome in c("Cases", "Hospitalization", "Death")) {
  for(vaccine_status in unique(waning_data$Vaccine.Status)) {
    for(prior_inf in unique(waning_data$Prior.Infection)) {
      subset_data <- waning_data %>% filter(Outcome == outcome,
                                            Vaccine.Status == vaccine_status)
      # ,
      #                                       Prior.Infection == prior_inf)
      
      if(length(unique(subset_data$Month3)) != length(subset_data$Month3)) {
        subset_data <- subset_data[1,]
        ve_data <- data.frame(ve = as.numeric(as.vector(subset_data[c(5:17)])),
                              months = c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
        
        
        p <- ggplot(ve_data, aes(x = months, y = ve)) +
          geom_line(size = 1.5) +
          ggtitle(paste0("Outcomes: ", outcome, "\nVaccine Status: ", vaccine_status, "\nPrior Infection: ", prior_inf, "\nNot Age-Stratified")) +
          ylab("Protective Effectiveness (%)")+
          ylim(0, 100) +
          scale_x_continuous(name="Months since vaccination/infections", limits=c(0, 12), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
        
       plot_list[[counter]] <- p
       counter <- counter + 1
      }
      
      else{
        melted_data <- melt(subset_data %>% select(Age, Month0.5:Month12)) %>%
          mutate(months = as.numeric(str_sub(as.character(variable), 6)))
        
        p <- ggplot(melted_data, aes(x = months, y = value, color = Age)) +
          geom_line(size = 1.5) +
          ggtitle(paste0("Outcomes: ", outcome, "\nVaccine Status: ", vaccine_status, "\nPrior Infection: ", prior_inf, "\nAge-Stratified")) +
          ylab("Protective Effectiveness (%)")+
          ylim(0, 100) +
          scale_x_continuous(name="Months since vaccination/infections", limits=c(0, 12), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
        plot_list[[counter]] <- p
        counter <- counter + 1
      }
      
    }
    
  }
  
}

plot_list = list()
counter <- 1
for (outcome in c("Cases", "Severe Outcomes")) {
  
  subset_data <- waning_data %>% filter(Outcome == outcome) %>% select(Vaccine.Status, Age, Prior.Infection, Month0.5:Month12)
  
  removing_duplicates <- subset_data %>% add_count(Month0.5, Month6, Month12) %>% distinct(Month0.5, Month6,Month12, .keep_all = TRUE) %>% filter(n < 6)
  
  melted_data <- merge(melt(removing_duplicates %>% select(-n), id_vars = c("Vaccine.Status", "Age", "Prior.Infection", "n")) %>%
                         mutate(months = as.numeric(str_sub(as.character(variable), 6))),
                       removing_duplicates %>% select(Vaccine.Status, Age, Prior.Infection, n)) %>%
    mutate(group = ifelse(n == 3, 
                          paste0(Vaccine.Status, ", ", Prior.Infection, " prior infection"),
                          paste0(Vaccine.Status, ", ", Age, ", ", Prior.Infection, " prior infection")))
  
  p <- ggplot(melted_data, aes(x = months, y = value, color = group)) +
    geom_line(size = 1) +
    ggtitle(paste0("Outcomes: ", outcome)) +
    ylab("Protective Effectiveness (%)")+
    ylim(0, 100) +
    scale_x_continuous(name="Months since vaccination/infections", limits=c(0, 12), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
    scale_color_manual(values = c(
      "dodgerblue2", "#E31A1C", # red
      "green4",
      "#6A3D9A", # purple
      "#FF7F00", # orange
      "black", "gold1",
      "skyblue2", "#FB9A99", # lt pink
      "palegreen2",
      "#CAB2D6", # lt purple
      "#FDBF6F", # lt orange
      "gray70", "khaki2",
      "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3",
      "darkorange4", "brown"
    ))
  
  plot_list[[counter]] <- p
  counter <- counter + 1

  
}



pdf("plots_combined_absoluteVE.pdf")
for (i in 1:(counter - 1)) {
  print(plot_list[[i]])
}
dev.off()



subset_data <- waning_data %>% filter(Outcome == "Severe Outcomes") %>% select(Vaccine.Status, Age, Prior.Infection, Month0.5:Month12)

removing_duplicates <- subset_data %>% add_count(Month0.5, Month3, Month6, Month12) %>% distinct(Month0.5, Month3,Month6,Month12, .keep_all = TRUE) %>% filter(n < 6,
                                                                                                                                                #case_when((Vaccine.Status %in% c("Unvaccinated", "2-dose")) | (Vaccine.Status == "3-dose" & Prior.Infection == "Yes") ~ Age == "18-49 years", 
                                                                                                                                                case_when((Vaccine.Status %in% c("Unvaccinated"))~ Age == "18-49 years",
                                                                                                                                                         T ~ Age == "18-49 years"))

melted_data <- merge(melt(removing_duplicates %>% select(-n), id_vars = c("Vaccine.Status", "Age", "Prior.Infection", "n")) %>%
                       mutate(months = as.numeric(str_sub(as.character(variable), 6))),
                     removing_duplicates %>% select(Vaccine.Status, Age, Prior.Infection, n)) %>%
  mutate(group = ifelse(n == 3, 
                        paste0(Vaccine.Status, ", ", Prior.Infection, " prior infection"),
                        paste0(Vaccine.Status, ", ", Age, ", ", Prior.Infection, " prior infection")))

ggplot(melted_data, aes(x = months, y = value, color = group)) +
  geom_line(size = 1) +
  ggtitle(paste0("Outcomes: ", "Severe Outcomes, Age: 18-49 years")) +
  ylab("Protective Effectiveness (%)")+
  ylim(0, 100) +
  scale_x_continuous(name="Months since vaccination/infections", limits=c(0, 12), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  scale_color_manual(values = c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  ))



#####################################
library(metafor)








