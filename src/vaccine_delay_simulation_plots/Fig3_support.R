prep_f3_data <- function(simulation_output, combined){
  sim_data <- simulation_output$state_by_age
  #There is no quantile data for age-strat
  sim_data <- filter(sim_data, region == 'england')
  sim_data <- subset(sim_data, select=-c(vaccine_daily_doses,
                                         booster_daily_doses,
                                         beta_step,
                                         scenario,
                                         analysis, quantile))
  sim_data <- filter(sim_data, state == 'diagnoses_admitted_inc')
  
  #Combine age and vaccine status:
  sim_data$category <- paste(sim_data$vaccine_status, sim_data$group, sep="_")
  
  #Load in fits age data
  fit_data <- combined$simulate$state_by_age$diagnoses_admitted
  fit_data <- aperm(fit_data, c(1,2,4,3))
  fit_data <- rowSums(fit_data, dims = 3L)
  for(i in 1:4){
    for(j in 1:4){
      fit_data[i,j,] <-c(0, diff(fit_data[i,j,]))
    }
  }
  fit_data <- fit_data[,,-1]
  
  age_classes <- rownames(fit_data[,,1])
  vacc_classes <- colnames(fit_data[,,1])
  
  fit_data_frame <- data.frame(state = NA, region = NA, date = as.Date("2020-01-01"),
                               value = NA, category = NA, group = NA, vaccine_status = NA)
  
  for(i in 1:length(age_classes)){
    for(j in 1:length(vacc_classes)){
      hold_data <- data.frame(state = rep('diagnoses_admitted_inc', 272),
                              region = rep('england', 272),
                              date = sircovid::sircovid_date_as_date(combined$simulate$date[2:273]),
                              value = fit_data[i,j,],
                              category = rep(paste(vacc_classes[j], age_classes[i], sep = '_'), 272),
                              group = rep(age_classes[i], 272),
                              vaccine_status = rep(vacc_classes[j], 272))
      
      fit_data_frame <- rbind(fit_data_frame, hold_data)
    }
  }
  
  #Remove the placeholder row
  fit_data_frame <- fit_data_frame[-1,]
  
  #We now want to make sure they start from 0.
  #Values are only initially present in the four unvaccinated categories
  sim_data <- filter(sim_data, date > as.Date("2020-12-15"))
  sim_data <- filter(sim_data, date < as.Date("2021-09-14"))
  
  fit_data_frame <- filter(fit_data_frame, date > as.Date("2020-12-15"))
  
  #Now combine the sim and fit data with separate column for each 
  names(sim_data)[names(sim_data) == 'value'] <- 'sim_value'
  names(fit_data_frame)[names(fit_data_frame) == 'value'] <- 'fit_value'
  
  var_names <- names(fit_data_frame)
  
  ret <- merge(sim_data, fit_data_frame, by = var_names[-4])
  
  ret$category <- as.factor(ret$category)
  
  
    ret <- pivot_longer(ret, c(sim_value, fit_value), names_to = "scenario", values_to = "value")
  
  return(ret)
}


plot_f3 <- function(f3_data){
  
  # #We prepare real-life data to fit to
  # admissions_deaths_vacc_status$date <- as.Date(admissions_deaths_vacc_status$date)
  # #First, we need to collapse age groups 
  # admissions_deaths_vacc_status$age_band_min <- as.factor(admissions_deaths_vacc_status$age_band_min)
  # init_levels <- levels(admissions_deaths_vacc_status$age_band_min)
  # new_levels <- c("age_0",  "age_0",  "age_0", "age_0", "age_0", "age_0", "age_0",
  #                  "age_35", "age_35", "age_35", "age_35",
  #                  "age_55", "age_55", "age_55", "age_55", "age_75", "age_75")
  # real_age_hosp_data <- admissions_deaths_vacc_status
  # real_age_hosp_data$age_band_min <- factor(new_levels[real_age_hosp_data$age_band_min])
  # real_age_hosp_data <- subset(real_age_hosp_data, select = -c(deaths))
  # real_age_hosp_data$vacc_status <- ifelse(is.na(real_age_hosp_data$vacc_status), "Unknown", real_age_hosp_data$vacc_status)
  # 
  # real_age_hosp_data <- aggregate(real_age_hosp_data$admissions,  list(date = real_age_hosp_data$date,
  #                                                                      age_band_min = real_age_hosp_data$age_band_min,
  #                                                                      vacc_status = real_age_hosp_data$vacc_status),  sum)
  # colnames(real_age_hosp_data) <- c("date", "age_band_min", "vacc_status", "admissions")
  # real_age_hosp_data <- aggregate(real_age_hosp_data$admissions,  list(date = real_age_hosp_data$date,
  #                                                                      age_band_min = real_age_hosp_data$age_band_min),  sum)
  # colnames(real_age_hosp_data) <- c("date", "group", "admissions")
  # real_age_hosp_data$scenario <- rep("fit_value", 1120)
  # 
  
  
  f3_data$vaccine_status <- factor(f3_data$vaccine_status, levels = c("unvaccinated", "waned_protection", "full_protection", "partial_protection" ) )
  
  type.labs <- c("3-week strategy", "12-week strategy")
  names(type.labs) <- c("sim_value", "fit_value")
  group.labs <- c("Age 0-34", "Age 35-54", "Age 55-74", "Age 75+")
  names(group.labs) <- c("age_0", "age_35", "age_55", "age_75")
  
  ggplot(data = f3_data,
         aes(x = date, y = value, fill = vaccine_status)) + 
    geom_area() + 
    facet_grid( rows = vars(scenario), cols = vars(group),
                              labeller = labeller(scenario = type.labs, group = group.labs)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y", name = 'Date') +
    labs(fill = "Vaccination status", y = "Daily hospital admissions", x = "Date") +
    scale_fill_manual(values = c("unvaccinated" ="lightsteelblue", "waned_protection" = "#26828E", "full_protection" = "#482878", "partial_protection" = "#FDE725"),
                      breaks = c("unvaccinated", "partial_protection", "full_protection", "waned_protection"),
                      labels = c(#"Unvaccinated", 
                                 #'One dose',
                                 "No vaccine \nprotection",
                                 "> 21 days post \nfirst dose", 
                                 'Full second dose \nprotection', 'Reduced second \ndose protection')
                      ) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(
            color="black", linetype="solid") 
            ) +
    #scale_fill_discrete(breaks=c('B', 'C', 'A')) +
    theme(legend.spacing.y = unit(0.65, 'cm'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE))
    
}






plot_f3_collapsed <- function(f3_data){
  
  #Collapse waned and 2 dose into same
  f3_data2 <- fct_collapse(f3_data$vaccine_status,
                           Unvaccinated = "unvaccinated",
                           One_Dose = "partial_protection",
                           Two_Doses = c("full_protection", "waned_protection"))
  
  f3_data$vaccine_status <- f3_data2
  
  #Aggregate over age and vaccine status
  
  aggregate(f3_data$value,  list(date = f3_data$date,
                                   group = f3_data$group,
                                   vaccine_status = f3_data$vaccine_status,
                                 scenario = f3_data$scenario
  ),  sum) -> f3_data
  
  colnames(f3_data) <- c("date", "group", "vaccine_status", "scenario", "value")
  
  type.labs <- c("3-week strategy", "12-week strategy")
  names(type.labs) <- c("sim_value", "fit_value")
  group.labs <- c("Age 0-34", "Age 35-54", "Age 55-74", "Age 75+")
  names(group.labs) <- c("age_0", "age_35", "age_55", "age_75")
  
  ggplot(data = f3_data,
         aes(x = date, y = value, fill = vaccine_status)) + 
    geom_area() + 
    facet_grid( rows = vars(scenario), cols = vars(group),
                labeller = labeller(scenario = type.labs, group = group.labs)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y", name = 'Date') +
    labs(fill = "Vaccination status", y = "Daily hospital admissions", x = "Date") +
    scale_fill_manual(values = c("Unvaccinated" ="lightsteelblue", "Two_Doses" = "#482878", "One_Dose" = "#FDE725"),
                      breaks = c("Unvaccinated", "One_Dose", "Two_Doses"),
                      #labels = c("Unvaccinated", 'One dose', 'Two doses')
                      labels = c("No vaccine \nprotection", '> 21 days post \nfirst dose',
                      '> 7 days post \nsecond dose')
    ) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(
            color="black", linetype="solid") 
    ) +
    #scale_fill_discrete(breaks=c('B', 'C', 'A')) +
    theme(legend.spacing.y = unit(0.65, 'cm'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE))
  
}



########################################
#SI Fits comparison


plot_f3_SI_stacked <- function(f3_data){
  ###Load in nation total data
  nation_2022 <- read.csv("nation_2022_05_03.csv")
  nation_2022$date <- as.Date(nation_2022$date)
  nation_2022 <- filter(nation_2022, date < as.Date("2021-09-14"))
  nation_2022 <- filter(nation_2022, date > as.Date("2020-12-14"))
  
  #Aggregate over age
  nation_2022 <- aggregate(nation_2022$value, list(date = nation_2022$date), sum)
  nation_2022$value <- c(0,(nation_2022$x[2:273] - nation_2022$x[1:272]))
  nation_2022 <- nation_2022[-1,-2]
  
  ####
  #Now load in NHS data
  nhs_data <- read.csv("estimated_admissions_by_age.csv")
  names(nhs_data)[names(nhs_data) == 'ï..date'] <- 'date'
  nhs_data <- nhs_data[,c(1,14,15:18)]
  colnames(nhs_data) <- c("date", "total", "age_0", 
                          "age_35", "age_55", "age_75")
  nhs_data$date <- as.Date(nhs_data$date, "%d/%m/%Y")
  nhs_data <- filter(nhs_data, date < as.Date("2021-09-14"))
  nhs_data <- filter(nhs_data, date > as.Date("2020-12-15"))
  
  #Calculate the necessary multipliers
  nhs_multipliers <- nation_2022$value/nhs_data$total
  for(i in 3:6){
    nhs_data[,i] <- round(nhs_data[,i]*nhs_multipliers)
  }
  
  #Separate age and boundary as a separate column
  nhs_data <- nhs_data[,-2]
  nhs_data <- pivot_longer(nhs_data,
                           c(age_0,
                             age_35,
                             age_55,
                             age_75),
                           names_to = "group", values_to = "admissions")
  
  
  nhs_data$scenario <- "fit_value"
  
  
  
  #Collapse f3 vaccine_status
  f3_data <- filter(f3_data, scenario == "fit_value")
  f3_data <- aggregate(f3_data$value, by = list(date = f3_data$date,
                                                group = f3_data$group), sum)
  colnames(f3_data) <- c("date", "group", "value")
  
 
  
  nhs_data$admissions[nhs_data$group == "age_55"] <- nhs_data$admissions[nhs_data$group == "age_55"] + nhs_data$admissions[nhs_data$group == "age_75"]
  nhs_data$admissions[nhs_data$group == "age_35"] <- nhs_data$admissions[nhs_data$group == "age_35"] + nhs_data$admissions[nhs_data$group == "age_55"]
  nhs_data$admissions[nhs_data$group == "age_0"] <- nhs_data$admissions[nhs_data$group == "age_0"] + nhs_data$admissions[nhs_data$group == "age_35"]
  
  #Convert to moving averages over 7 days
  moving_average <- function(x, n = 7) {             
    stats::filter(x, rep(1 / n, n), sides = 2)
  }
  
  nhs_data$admissions[nhs_data$group == "age_0"] <- moving_average(nhs_data$admissions[nhs_data$group == "age_0"])
  nhs_data$admissions[nhs_data$group == "age_35"] <- moving_average(nhs_data$admissions[nhs_data$group == "age_35"])
  nhs_data$admissions[nhs_data$group == "age_55"] <- moving_average(nhs_data$admissions[nhs_data$group == "age_55"])
  nhs_data$admissions[nhs_data$group == "age_75"] <- moving_average(nhs_data$admissions[nhs_data$group == "age_75"])
  
  
  ggplot(data = f3_data,
         aes(x = date, fill = group)) + 
    geom_area(aes(y = value), alpha = 0.8) + 
    geom_line(data = nhs_data, aes(x = date, y = admissions,
                                                color = group),
               alpha = 1, size = 1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y", name = 'Date', expand = c(0, 0)) +
    labs(fill = "Model output \nage groups", y = "Daily hospital admissions", x = "Date", color = "NHS data \nage groups") +
    scale_fill_manual(values = colorspace::lighten(c("#841744", "#26828E", "#482878", "#FDE725"),0.2),
                      labels = c("0-34", '34-54', '55-74', '75+')) +
    scale_color_manual(values = c( colorspace::darken("#841744", 0.4),
                                   colorspace::darken("#26828E",0.4),
                                   colorspace::darken("#482878",0.4),
                                   colorspace::darken("#FDE725",0.7) ),
                      labels = c("0-34", '35-54', '55-74', '75+')) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="solid")
    ) 
  
}



plot_f3_SI_age_prop <- function(f3_data){
  
  f3_data$vaccine_status <- factor(f3_data$vaccine_status, levels = c("unvaccinated", "waned_protection", "full_protection", "partial_protection" ) )
  
  type.labs <- c("3-week strategy", "12-week strategy")
  names(type.labs) <- c("sim_value", "fit_value")
  group.labs <- c("0-34", "35-54", "55-74", "75+")
  names(group.labs) <- c("age_0", "age_35", "age_55", "age_75")
  
  f3_data_test <- f3_data
  
  f3_data_test <- aggregate(f3_data_test$value, by = list(date = f3_data_test$date,
                                                          group = f3_data_test$group,
                                                          scenario = f3_data_test$scenario), sum)
  colnames(f3_data_test) <- c("date", "group", "scenario", "value")
  
  f3_data_test2 <- f3_data_test
  
  Dates_list <- unique(f3_data_test$date)
  group_list <- unique(f3_data_test$group)
  for(I in 1:length(Dates_list)){
    i <- Dates_list[I]
    for(J in 1:length(group_list)){
      j <- group_list[J]
      f3_data_test2$value[f3_data_test$date == i & f3_data_test$group == j & f3_data_test$scenario == "fit_value"] <- f3_data_test$value[f3_data_test$date == i & f3_data_test$group == j & f3_data_test$scenario == "fit_value"] / sum(f3_data_test$value[f3_data_test$date == i & f3_data_test$scenario == "fit_value"])
      f3_data_test2$value[f3_data_test$date == i & f3_data_test$group == j & f3_data_test$scenario == "sim_value"] <- f3_data_test$value[f3_data_test$date == i & f3_data_test$group == j & f3_data_test$scenario == "sim_value"] / sum(f3_data_test$value[f3_data_test$date == i & f3_data_test$scenario == "sim_value"])
      
    }
    
  }
  
  
  
  ggplot(data = f3_data_test2,
         aes(x = date, y = value, fill = group)) + 
    geom_area() + facet_grid( rows = vars(scenario),
                              labeller = labeller(scenario = type.labs)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y", name = 'Date') +
    labs(fill = "Age group", y = "Proportion of hospital admissions", x = "Date") +
    scale_fill_manual(values = colorspace::lighten(c("#841744", "#26828E", "#482878", "#FDE725"),0.2),
                      labels = c("0-34", '35-54', '55-74', '75+')) +
    theme(axis.text=element_text(size=rel(1.2)),
          panel.spacing = unit(1, "lines"),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.background = element_rect(
            color="black", linetype="solid") 
    )
  
}


####################################################
#  Hospitalisations presented by age breakdown in model fits

plot_f3_SI_hosp_dist <- function(f3_data){
  
  hold_data <- f3_data
  #Collapse waned and 2 dose into same
  f3_data2 <- fct_collapse(f3_data$vaccine_status,
                          Unvaccinated = "unvaccinated",
                          One_Dose = "partial_protection",
                          Two_Doses = c("full_protection", "waned_protection"))
  
  hold_data$vaccine_status <- f3_data2
  
  #Remove sim values
  hold_data <- filter(hold_data, scenario == "fit_value")
  #Aggregate over age and vaccine status
  
  aggregate(hold_data$value,  list(date = hold_data$date,
                                 group = hold_data$group,
                                 vaccine_status = hold_data$vaccine_status
                                 ),  sum) -> hold_data
  aggregate(hold_data$x, list(date = hold_data$date,
                              vaccine_status = hold_data$vaccine_status
  ),  sum) -> hold_data
  
  colnames(hold_data) <- c("date", "vaccine_status", "value")
  
 #Collect totals to calculate percentage from
  total_hold <- aggregate(hold_data$value,
                          list(data = hold_data$date),
                          sum)
  colnames(total_hold) <- c("date", "percentage")
  
  hold_data <- left_join(hold_data, total_hold)
  
  hold_data$percentage <- hold_data$value/hold_data$percentage
  
  
  ggplot(data = hold_data,
         aes(x = date, y = value, fill = vaccine_status)) + 
    geom_area() + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y", name = 'Date') +
    labs(fill = "Vaccination status", y = "Daily hospital admissions", x = "Date") +
    scale_fill_manual(values = c("Unvaccinated" ="lightsteelblue", "Two_Doses" = "#482878", "One_Dose" = "#FDE725"),
                      breaks = c("Unvaccinated", "One_Dose", "Two_Doses"),
                      #labels = c("Unvaccinated", 'One dose', 'Two doses')
                      labels = c("No vaccine \nprotection", '> 21 days post \nfirst dose',
                                 '> 7 days post \nsecond dose')
    ) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)) 
    ) +
    #scale_fill_discrete(breaks=c('B', 'C', 'A')) +
    theme(legend.spacing.y = unit(0.65, 'cm'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE)) -> p1
  
  ggplot(data = hold_data,
         aes(x = date, y = percentage, color = vaccine_status)) + 
    geom_line() + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b %y", name = 'Date') +
    labs(color = "Vaccination status", y = "Daily hospital admissions", x = "Date") +
    scale_color_manual(values = c("Unvaccinated" ="lightsteelblue", "Two_Doses" = "#482878", "One_Dose" = "#FDE725"),
                      breaks = c("Unvaccinated", "One_Dose", "Two_Doses"),
                      labels = c("Unvaccinated", 'One dose', 'Two doses')
    ) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)) 
    ) +
    #scale_fill_discrete(breaks=c('B', 'C', 'A')) +
    theme(legend.spacing.y = unit(0.65, 'cm'))  +
    ## important additional element
    guides(color = guide_legend(byrow = TRUE)) -> p2
  
  plot_grid(p1,p2, nrow = 2)
  
}




