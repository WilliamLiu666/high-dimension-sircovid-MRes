prep_f2p1_data <- function(simulation_output, combined){
  
  #Extract n_doses for each
  sim_n_doses <- simulation_output$n_doses
  fit_n_doses_raw <- combined$simulate$n_doses
  fit_n_doses_raw <- fit_n_doses_raw[,,8,]
  
  sim_n_doses <- filter(sim_n_doses, region == "england" & state != "booster_dose")
  sim_n_doses <- filter(sim_n_doses, date < as.Date("2021-09-14"))
  sim_n_doses <- subset(sim_n_doses, select=-c(vaccine_daily_doses,
                                                    booster_daily_doses,
                                                    beta_step,
                                                    scenario, quantile,
                                                    analysis, vaccine_status))

  n_doses_states <- dimnames(fit_n_doses_raw)[[2]]
  n_doses_ages <- unique(sim_n_doses$group)
  fit_n_doses <- data.frame(state = NA, region = NA, group = NA,
                            date = as.Date("2000-02-02") #placeholder!
                            , value = NA)
  #Age
  for(i in 1:19){
    #state
    for(j in 1:8){
      hold <- data.frame(state = rep(n_doses_states[j], 273), 
                         region = rep("england", 273), 
                         group = rep(n_doses_ages[i], 273),
                         date = sircovid::sircovid_date_as_date(combined$simulate$date), 
                         value = fit_n_doses_raw[i,j,])
      
      fit_n_doses <- rbind(fit_n_doses, hold)
    }
  }
  
  sim_n_doses$type <- rep("counterfactual", length(sim_n_doses$state))
  fit_n_doses$type <- rep("fit", length(fit_n_doses$state))
  #Remove the placeholder row
  fit_n_doses <- fit_n_doses[-1,]
  
  All_n_doses <- rbind(sim_n_doses, fit_n_doses)
  
  
  #Lastly, we aggregate by age
  All_n_doses <- aggregate(All_n_doses$value, list(state = All_n_doses$state,
                                             region = All_n_doses$region,
                                             date = All_n_doses$date,
                                             type = All_n_doses$type), sum)
  colnames(All_n_doses) <- c("state", "region", "date", "type", "value")
  
  return(All_n_doses)
  
}

#############
#PLOTTING
#############

plot_f2p1 <- function(f2p1_data){
  
  #First, we wish to only look at cumulative incidence
  keep_states <- c("first_dose", "second_dose", "waned_dose")
  
  f2p1_data <- filter(f2p1_data, state %in% keep_states)  
  
  #We will need to do some recalculating here too. As 1st dose cumulative should have
  #2nd dose removed from it, and waned dose should then be removed from second dose
  f2p1_data <- filter(f2p1_data, date > as.Date("2020-12-14"))
  f2p1_data_proc <- f2p1_data
  
  for(i in 1:length(f2p1_data_proc$state)){
    if(f2p1_data_proc$date[i] != as.Date("2020-12-15")){
      
      hold_state <- f2p1_data_proc$state[i]
      hold_date <- f2p1_data_proc$date[i]
      hold_type <- f2p1_data_proc$type[i]
      
      if(hold_state == "first_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "second_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }else if(hold_state == "second_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "waned_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }
      
    }
  }
  
  #We now want to add in a new category, unvaccinated that takes the total up to population:
  england_pop <- sum(sircovid:::sircovid_population("england"))
  f2p1_data_proc_new <- f2p1_data_proc
  dates <- as.Date(unique(f2p1_data_proc$date))
  for(I in 1:length(dates)){
    i <- dates[I]
    for(j in c("fit", "counterfactual")){
      hold_data <- data.frame(state = "unvaccinated", 
                              region = "england",
                              date = i,
                              type = j,
                              value = england_pop -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "first_dose"])
                              -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "second_dose"]) 
                              -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "waned_dose"])
      )
      f2p1_data_proc_new <- rbind(f2p1_data_proc_new, hold_data)
      
    }
  }
  
  f2p1_data_proc <- f2p1_data_proc_new
  type.labs <- c("Counterfactual", "Model fit")
  names(type.labs) <- c("counterfactual", "fit")
  f2p1_data_proc$type = factor(f2p1_data_proc$type, levels=c("fit", "counterfactual"))
  
  f2p1_data_proc$state <- factor(f2p1_data_proc$state, levels = c("unvaccinated", "first_dose",
                                                                  "second_dose",
                                                                  "waned_dose"
                                                                  ) #rev(levels(f2p1_data_proc$state))
                                 )
  
  #Filter down to just first of each month
  keep_dates <- seq(as.Date("2021/2/1"), by = "month", length.out = 9)
  
  f2p1_data_proc_new <- filter(f2p1_data_proc, date %in% keep_dates)
  
  
  ggplot(data = f2p1_data_proc_new, aes(x = date, y = value/1000000, fill = state)) +
    geom_col(data = filter(f2p1_data_proc_new, type == "fit"), aes(x = date-6,group = type), width = 10) +
    geom_text(data = filter(f2p1_data_proc_new, type == "fit" & state == "first_dose"), aes(label=rep("12w", 8),x = date-6, y = 58), color = "#B25A82", fontface = "bold") +
    geom_col(data = filter(f2p1_data_proc_new, type == "counterfactual"), aes(x = date+6, group = type), width = 10) +
    geom_text(data = filter(f2p1_data_proc_new, type == "counterfactual" & state == "first_dose"), aes(label=rep("3w", 8),x = date+6, y = 58), color = "#9BC362", fontface = "bold") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,59)) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
    labs(fill = "Vaccination status", y = "Population \n(Millions)", x = "") +
    #VIRIDIS
    #"lavenderblush2"
    scale_fill_manual(values = c("lightsteelblue", "#FDE725", "#482878", "#26828E"),
    
                      #ARCHAMBAULT
    #scale_fill_manual(values = c("gray80", "#88a0dc", "#7c4b73", "#ed968c"),
    #MUTED
    #scale_fill_manual(values = c("gray80", "#CC6677", "#332288", "#DDCC77" ),
    #KANDINSKY
    #scale_fill_manual(values = rev(c("#3b7c70", "#ce9642", "#898e9f", "#3b3a3e")),
    
                    labels = c("Unvaccinated", 'One Dose',
                               'Full second dose \nprotection',
                               'Reduced second \ndose protection'
                                  )) +
   
   theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="blank")
    ) +
    theme(legend.spacing.y = unit(0.4, 'cm'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE))-> p1
  
  

  
  return(p1)
  
}




plot_f2p1_collapsed <- function(f2p1_data){
  
  #First, we wish to only look at cumulative incidence
  keep_states <- c("first_dose", "second_dose", "waned_dose")
  
  f2p1_data <- filter(f2p1_data, state %in% keep_states)  
  
  #We will need to do some recalculating here too. As 1st dose cumulative should have
  #2nd dose removed from it, and waned dose should then be removed from second dose
  f2p1_data <- filter(f2p1_data, date > as.Date("2020-12-14"))
  f2p1_data_proc <- f2p1_data
  
  for(i in 1:length(f2p1_data_proc$state)){
    if(f2p1_data_proc$date[i] != as.Date("2020-12-15")){
      
      hold_state <- f2p1_data_proc$state[i]
      hold_date <- f2p1_data_proc$date[i]
      hold_type <- f2p1_data_proc$type[i]
      
      if(hold_state == "first_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "second_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }else if(hold_state == "second_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "waned_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }
      
    }
  }
  
  #We now want to add in a new category, unvaccinated that takes the total up to population:
  england_pop <- sum(sircovid:::sircovid_population("england"))
  f2p1_data_proc_new <- f2p1_data_proc
  dates <- as.Date(unique(f2p1_data_proc$date))
  for(I in 1:length(dates)){
    i <- dates[I]
    for(j in c("fit", "counterfactual")){
      hold_data <- data.frame(state = "unvaccinated", 
                              region = "england",
                              date = i,
                              type = j,
                              value = england_pop -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "first_dose"])
                              -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "second_dose"]) 
                              -(f2p1_data_proc$value[f2p1_data_proc$date == i & f2p1_data_proc$type == j & f2p1_data_proc$state == "waned_dose"])
      )
      f2p1_data_proc_new <- rbind(f2p1_data_proc_new, hold_data)
      
    }
  }
  
  f2p1_data_proc <- f2p1_data_proc_new
  type.labs <- c("Counterfactual", "Model fit")
  names(type.labs) <- c("counterfactual", "fit")
  f2p1_data_proc$type = factor(f2p1_data_proc$type, levels=c("fit", "counterfactual"))
  
  f2p1_data_proc$state <- factor(f2p1_data_proc$state, levels = c("unvaccinated", "first_dose",
                                                                  "second_dose",
                                                                  "waned_dose"
  ) #rev(levels(f2p1_data_proc$state))
  )
  
  #Filter down to just first of each month
  keep_dates <- seq(as.Date("2021/2/1"), by = "month", length.out = 9)
  
  f2p1_data_proc_new <- filter(f2p1_data_proc, date %in% keep_dates)
  
  
  #Collapse waned and 2 dose into same
  f2p1_data_proc_new2 <- fct_collapse(f2p1_data_proc_new$state,
                           Unvaccinated = "unvaccinated",
                           One_Dose = "first_dose",
                           Two_Doses = c("second_dose", "waned_dose"))
  
  f2p1_data_proc_new$state <-f2p1_data_proc_new2
  
  #Aggregate over age and vaccine status
  
  aggregate(f2p1_data_proc_new$value,  list(date = f2p1_data_proc_new$date,
                                 state = f2p1_data_proc_new$state,
                                 type = f2p1_data_proc_new$type
  ),  sum) -> f2p1_data_proc_new
  
  colnames(f2p1_data_proc_new) <- c("date", "state", "type","value")
  
  f2p1_data_proc_new <- f2p1_data_proc_new[nrow(f2p1_data_proc_new):1,]
  
  ggplot(data = f2p1_data_proc_new, aes(x = date, y = value/1000000, fill = state)) +
    geom_col(data = filter(f2p1_data_proc_new, type == "fit"),
             aes(x = date-6,group = type), width = 10) +
    geom_text(data = filter(f2p1_data_proc_new, type == "fit" & state == "One_Dose"), aes(label=rep("12w", 8),x = date-6, y = 58), color = "#B25A82", fontface = "bold") +
    geom_col(data = filter(f2p1_data_proc_new, type == "counterfactual"),
             aes(x = date+6, group = type), width = 10) +
    geom_text(data = filter(f2p1_data_proc_new, type == "counterfactual" & state == "One_Dose"), aes(label=rep("3w", 8),x = date+6, y = 58), color = "#9BC362", fontface = "bold") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,59)) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
    labs(fill = "Vaccination status", y = "Population \n(Millions)", x = "") +
    #VIRIDIS
    
    scale_fill_manual(values = c("lightsteelblue", #"lavenderblush2"
                                 "#FDE725", "#482878"),
                      
                      #ARCHAMBAULT
                      #scale_fill_manual(values = c("gray80", "#88a0dc", "#7c4b73", "#ed968c"),
                      #MUTED
                      #scale_fill_manual(values = c("gray80", "#CC6677", "#332288", "#DDCC77" ),
                      #KANDINSKY
                      #scale_fill_manual(values = rev(c("#3b7c70", "#ce9642", "#898e9f", "#3b3a3e")),
                      
                      #labels = c("Unvaccinated", 'One dose',
                      #           'Two doses'
                      #Technically we want:
                      labels = c("No vaccine \nprotection", '> 21 days post \nfirst dose',
                                            '> 7 days post \nsecond dose'
                      )) +
    
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="blank")
    ) +
    theme(legend.spacing.y = unit(0.4, 'cm'))  +
    ## important additional element
    guides(fill = guide_legend(byrow = TRUE))-> p1
  
  
  
  
  return(p1)
  
}

plot_f2p2 <- function(f2p2_data){
  
  #This plot will demonstrate population average VE over time.
  
  
  #First load the VE values
  #And average them across all of england
  regions <- sircovid::regions("england")
  average_vacc_efficacy <- combined$base$east_of_england$average_vacc_efficacy_alpha$central
  average_vacc_efficacy$death <- 0
  average_vacc_efficacy$severe_disease <- 0
  average_vacc_efficacy$infection <- 0
  for(i in 1:7){
    hold_efficacy <- combined$base[[i]]$average_vacc_efficacy_alpha$central
    average_vacc_efficacy$death <- average_vacc_efficacy$death + hold_efficacy$death*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$severe_disease <- average_vacc_efficacy$severe_disease + hold_efficacy$severe_disease*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$infection <- average_vacc_efficacy$infection + hold_efficacy$infection*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
  }
  
  #Now we want to average over all age groups for VE
  #Lose stratum5 and 1
  average_vacc_efficacy <- subset(average_vacc_efficacy, select = -c(analysis,
                                                                     prop_pfizer,
                                                                     disease,
                                                                     transmission))
  average_vacc_efficacy <- filter(average_vacc_efficacy, stratum %in% c("stratum_2", "stratum_3", "stratum_4"))
  #We also need to lose 18/19 group (CHW/CHR), as this population number is not replenished over time.
  average_vacc_efficacy <- filter(average_vacc_efficacy, group %in% 1:17)
  stratums <- unique(average_vacc_efficacy$stratum)
  average_vacc_efficacy_alpha <- data.frame(stratum = NA, 
                                            death = NA,
                                            severe_disease = NA,
                                            infection = NA)
  for(i in 1:length(stratums)){
    I <- stratums[i]
    hold_efficacy <- filter(average_vacc_efficacy, stratum == I)
    hold_efficacy$death <- hold_efficacy$death *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$severe_disease <- hold_efficacy$severe_disease *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$infection <- hold_efficacy$infection *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    
    hold <- data.frame(stratum = I,
                       death = sum(hold_efficacy$death),
                       severe_disease = sum(hold_efficacy$severe_disease),
                       infection = sum(hold_efficacy$infection))
    
    average_vacc_efficacy_alpha <- rbind(average_vacc_efficacy_alpha, hold)
    
  }
  average_vacc_efficacy_alpha <- average_vacc_efficacy_alpha[-1,]
  
  #DELTA TOO
  ############################
  average_vacc_efficacy <- combined$base$east_of_england$average_vacc_efficacy_delta$central
  average_vacc_efficacy$death <- 0
  average_vacc_efficacy$severe_disease <- 0
  average_vacc_efficacy$infection <- 0
  for(i in 1:7){
    hold_efficacy <- combined$base[[i]]$average_vacc_efficacy_delta$central
    average_vacc_efficacy$death <- average_vacc_efficacy$death + hold_efficacy$death*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$severe_disease <- average_vacc_efficacy$severe_disease + hold_efficacy$severe_disease*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$infection <- average_vacc_efficacy$infection + hold_efficacy$infection*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
  }
  
  #Now we want to average over all age groups for VE
  #Lose stratum5 and 1
  average_vacc_efficacy <- subset(average_vacc_efficacy, select = -c(analysis,
                                                                     prop_pfizer,
                                                                     disease,
                                                                     transmission))
  average_vacc_efficacy <- filter(average_vacc_efficacy, stratum %in% c("stratum_2", "stratum_3", "stratum_4"))
  #We also need to lose 18/19 group (CHW/CHR),
  average_vacc_efficacy <- filter(average_vacc_efficacy, group %in% 1:17)
  stratums <- unique(average_vacc_efficacy$stratum)
  average_vacc_efficacy_delta <- data.frame(stratum = NA, 
                                            death = NA,
                                            severe_disease = NA,
                                            infection = NA)
  for(i in 1:length(stratums)){
    I <- stratums[i]
    hold_efficacy <- filter(average_vacc_efficacy, stratum == I)
    hold_efficacy$death <- hold_efficacy$death *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$severe_disease <- hold_efficacy$severe_disease *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$infection <- hold_efficacy$infection *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    
    hold <- data.frame(stratum = I,
                       death = sum(hold_efficacy$death),
                       severe_disease = sum(hold_efficacy$severe_disease),
                       infection = sum(hold_efficacy$infection))
    
    average_vacc_efficacy_delta <- rbind(average_vacc_efficacy_delta, hold)
    
  }
  average_vacc_efficacy_delta <- average_vacc_efficacy_delta[-1,]
  
  
  #We now have all the VE, now we need to build the data frame we will plot
  #There will be six categories: 3*2 = (infection/hosp/death)*(fit/counterfactual)
  
  ###
  #Start by filtering out data to be as used in f2p1
  
  #First, we wish to only look at cumulative incidence
  keep_states <- c("first_dose", "second_dose", "waned_dose")
  
  f2p2_data <- filter(f2p2_data, state %in% keep_states)  
  
  #We will need to do some recalculating here too. As 1st dose cumulative should have
  #2nd dose removed from it, and waned dose should then be removed from second dose
  f2p2_data <- filter(f2p2_data, date > as.Date("2020-12-14"))
  f2p1_data_proc <- f2p2_data
  
  for(i in 1:length(f2p1_data_proc$state)){
    if(f2p1_data_proc$date[i] != as.Date("2020-12-15")){
      
      hold_state <- f2p1_data_proc$state[i]
      hold_date <- f2p1_data_proc$date[i]
      hold_type <- f2p1_data_proc$type[i]
      
      if(hold_state == "first_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "second_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }else if(hold_state == "second_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "waned_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }
      
    }
  }
  
  f2p1_data_proc$value <-  f2p1_data_proc$value/sum(sircovid:::sircovid_population("england"))
  #Now we make the new data frame that will be used for the fig, VE over time
  
  #JUST ALPHA HERE
  Alpha_VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "fit")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_alpha[,j])*(hold_population$value)))
      
      Alpha_VE_time_data_fit <- rbind(Alpha_VE_time_data_fit, hold)
    }
  }
  Alpha_VE_time_data_fit <- Alpha_VE_time_data_fit[-1,]
  
  
  #SAME BUT FOR COUNTERFACTUAL
  Alpha_VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "counterfactual")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_alpha[,j])*(hold_population$value)))
      
      Alpha_VE_time_data_counterfactual <- rbind(Alpha_VE_time_data_counterfactual, hold)
    }
  }
  Alpha_VE_time_data_counterfactual <- Alpha_VE_time_data_counterfactual[-1,]
  
  
  #Now for Delta
  Delta_VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "fit")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_delta[,j])*(hold_population$value)))
      
      Delta_VE_time_data_fit <- rbind(Delta_VE_time_data_fit, hold)
    }
  }
  Delta_VE_time_data_fit <- Delta_VE_time_data_fit[-1,]
  
  
  #SAME BUT FOR COUNTERFACTUAL
  Delta_VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "counterfactual")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_delta[,j])*(hold_population$value)))
      
      Delta_VE_time_data_counterfactual <- rbind(Delta_VE_time_data_counterfactual, hold)
    }
  }
  Delta_VE_time_data_counterfactual <- Delta_VE_time_data_counterfactual[-1,]
  
  
  
  
  
  #Stick them together
  Delta_VE_time_data_fit$type <- "fit"
  Delta_VE_time_data_counterfactual$type <- "counterfactual"
  Alpha_VE_time_data_fit$type <- "fit"
  Alpha_VE_time_data_counterfactual$type <- "counterfactual"
  
  #Load in the Delta emergence date
  Delta_prop <- read.csv("Delta_emergence_data.csv")
  Delta_prop <- aggregate(Delta_prop$pos_mean, by = list(dates = Delta_prop$dates), FUN = mean)
  
  ####
  #And can directly extract the same for the counterfactual (in case it differs)
  
  n_strains <- filter(simulation_output$state, state %in% c("n_strain_1", "n_strain_2"))
  n_strains <- filter(n_strains, region == "england")
  n_strains <- filter(n_strains, quantile == "50%")
  n_strain_1 <- filter(n_strains, state == "n_strain_1")
  n_strain_2 <- filter(n_strains, state == "n_strain_2")
  
  n_strain_1$value <- c(0, diff(n_strain_1$value))
  n_strain_2$value <- c(0, diff(n_strain_2$value))
  n_strain_1 <- n_strain_1[-1,]
  n_strain_2 <- n_strain_2[-1,]
  
  n_strains <- rbind(n_strain_1, n_strain_2)
  
  simulation_delta_proportion <- data.frame(date = n_strain_1$date,
                                            prop_delta = n_strain_2$value/(n_strain_1$value + n_strain_2$value))
  
  simulation_delta_proportion <- filter(simulation_delta_proportion, date > as.Date("2021-03-07"))
  simulation_delta_proportion <- filter(simulation_delta_proportion, date < as.Date("2021-08-01"))
  ####
  
  #Build the final data frame
  VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA, type = "fit")
  for(i in 1:length(dates)){
    I <- dates[i]
    if(I < as.Date("2021-03-08")){
      hold_Alpha <- filter(Alpha_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Alpha$VE[1]),
                                (hold_Alpha$VE[2]),
                                (hold_Alpha$VE[3])),
                         type = rep("fit", 3))
    } else if(I < as.Date("2021-08-01")){
      
      hold_Delta_prop <- Delta_prop$x[which(Delta_prop$dates == I)]/100
      hold_Alpha <- filter(Alpha_VE_time_data_fit, date == I )
      hold_Delta <- filter(Delta_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]*(hold_Delta_prop)) + (hold_Alpha$VE[1]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[2]*(hold_Delta_prop)) + (hold_Alpha$VE[2]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[3]*(hold_Delta_prop)) + (hold_Alpha$VE[3]*(1-hold_Delta_prop))),
                         type = rep("fit", 3))
      
    }else{
      
      hold_Delta <- filter(Delta_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]),
                                (hold_Delta$VE[2]),
                                (hold_Delta$VE[3])),
                         type = rep("fit", 3))
    }
    VE_time_data_fit <- rbind(VE_time_data_fit, hold)
    
  }
  VE_time_data_fit <- VE_time_data_fit[-1,]
  
  #REPEAT FOR COUNTERFACTUAL
  VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA, type = "counterfactual")
  for(i in 1:length(dates)){
    I <- dates[i]
    if(I < as.Date("2021-03-08")){
      hold_Alpha <- filter(Alpha_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Alpha$VE[1]),
                                (hold_Alpha$VE[2]),
                                (hold_Alpha$VE[3])),
                         type = rep("counterfactual", 3))
    } else if(I < as.Date("2021-08-01")){
      
      hold_Delta_prop <- simulation_delta_proportion$prop_delta[which(simulation_delta_proportion$date == I)]
      hold_Alpha <- filter(Alpha_VE_time_data_counterfactual, date == I )
      hold_Delta <- filter(Delta_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]*(hold_Delta_prop)) + (hold_Alpha$VE[1]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[2]*(hold_Delta_prop)) + (hold_Alpha$VE[2]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[3]*(hold_Delta_prop)) + (hold_Alpha$VE[3]*(1-hold_Delta_prop))),
                         type = rep("counterfactual", 3))
      
    }else{
      
      hold_Delta <- filter(Delta_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]),
                                (hold_Delta$VE[2]),
                                (hold_Delta$VE[3])),
                         type = rep("counterfactual", 3))
    }
    VE_time_data_counterfactual <- rbind(VE_time_data_counterfactual, hold)
    
  }
  VE_time_data_counterfactual <- VE_time_data_counterfactual[-1,]
  
  VE_time_data <- rbind(VE_time_data_fit, VE_time_data_counterfactual)
  
  #NOW we filter down for the stacked bars we want
  keep_dates <- seq(as.Date("2021/2/1"), by = "month", length.out = 9)
  VE_time_data_new <- filter(VE_time_data, date %in% keep_dates)
  
  #We can drop death
  VE_time_data_new <- filter(VE_time_data_new, outcome != "death")
  
  outcome.labs <- c("Infection", "Severe disease")
  names(outcome.labs) <- c("infection", "severe_disease")
  
  VE_time_data_new <- filter(VE_time_data_new, outcome == 'infection')
  ggplot(data = VE_time_data_new,
         aes(x = date, y = (1-VE))) +
    geom_col(data = filter(VE_time_data_new, type == "fit"), aes(x = date-6,group = type, fill = type), width = 10) +
    geom_col(data = filter(VE_time_data_new, type == "counterfactual"), aes(x = date+6, group = type, fill = type), width = 10) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0.4,1), oob = rescale_none) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
    labs(fill = "Strategy", y = "Proportion at risk of infection", x = "") +
    geom_hline(yintercept = c(0.6, 0.8, 1), color = "black", alpha = 0.15) +
    #facet_grid(rows = vars(outcome),
    #           labeller = labeller(outcome = outcome.labs)) +
    scale_fill_manual(breaks = c("fit", "counterfactual"),
                      values = c("#B25A82", "#9BC362"),
                      labels = c("12-week", '3-week')) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="blank")
    ) 
  
  
}






plot_f2p2_inf_and_hosp <- function(f2p2_data){
  
  #This plot will demonstrate population average VE over time.
  
  
  #First load the VE values
  #And average them across all of england
  regions <- sircovid::regions("england")
  average_vacc_efficacy <- combined$base$east_of_england$average_vacc_efficacy_alpha$central
  average_vacc_efficacy$death <- 0
  average_vacc_efficacy$severe_disease <- 0
  average_vacc_efficacy$infection <- 0
  for(i in 1:7){
    hold_efficacy <- combined$base[[i]]$average_vacc_efficacy_alpha$central
    average_vacc_efficacy$death <- average_vacc_efficacy$death + hold_efficacy$death*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$severe_disease <- average_vacc_efficacy$severe_disease + hold_efficacy$severe_disease*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$infection <- average_vacc_efficacy$infection + hold_efficacy$infection*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
  }
  
  #Now we want to average over all age groups for VE
  #Lose stratum5 and 1
  average_vacc_efficacy <- subset(average_vacc_efficacy, select = -c(analysis,
                                                                     prop_pfizer,
                                                                     disease,
                                                                     transmission))
  average_vacc_efficacy <- filter(average_vacc_efficacy, stratum %in% c("stratum_2", "stratum_3", "stratum_4"))
  #We also need to lose 18/19 group (CHW/CHR), as this population number is not replenished over time.
  average_vacc_efficacy <- filter(average_vacc_efficacy, group %in% 1:17)
  stratums <- unique(average_vacc_efficacy$stratum)
  average_vacc_efficacy_alpha <- data.frame(stratum = NA, 
                                            death = NA,
                                            severe_disease = NA,
                                            infection = NA)
  for(i in 1:length(stratums)){
    I <- stratums[i]
    hold_efficacy <- filter(average_vacc_efficacy, stratum == I)
    hold_efficacy$death <- hold_efficacy$death *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$severe_disease <- hold_efficacy$severe_disease *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$infection <- hold_efficacy$infection *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    
    hold <- data.frame(stratum = I,
                       death = sum(hold_efficacy$death),
                       severe_disease = sum(hold_efficacy$severe_disease),
                       infection = sum(hold_efficacy$infection))
    
    average_vacc_efficacy_alpha <- rbind(average_vacc_efficacy_alpha, hold)
    
  }
  average_vacc_efficacy_alpha <- average_vacc_efficacy_alpha[-1,]
  
  #DELTA TOO
  ############################
  average_vacc_efficacy <- combined$base$east_of_england$average_vacc_efficacy_delta$central
  average_vacc_efficacy$death <- 0
  average_vacc_efficacy$severe_disease <- 0
  average_vacc_efficacy$infection <- 0
  for(i in 1:7){
    hold_efficacy <- combined$base[[i]]$average_vacc_efficacy_delta$central
    average_vacc_efficacy$death <- average_vacc_efficacy$death + hold_efficacy$death*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$severe_disease <- average_vacc_efficacy$severe_disease + hold_efficacy$severe_disease*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
    average_vacc_efficacy$infection <- average_vacc_efficacy$infection + hold_efficacy$infection*(sum(sircovid:::sircovid_population(regions[i]))/sum(sircovid:::sircovid_population("england")))
  }
  
  #Now we want to average over all age groups for VE
  #Lose stratum5 and 1
  average_vacc_efficacy <- subset(average_vacc_efficacy, select = -c(analysis,
                                                                     prop_pfizer,
                                                                     disease,
                                                                     transmission))
  average_vacc_efficacy <- filter(average_vacc_efficacy, stratum %in% c("stratum_2", "stratum_3", "stratum_4"))
  #We also need to lose 18/19 group (CHW/CHR),
  average_vacc_efficacy <- filter(average_vacc_efficacy, group %in% 1:17)
  stratums <- unique(average_vacc_efficacy$stratum)
  average_vacc_efficacy_delta <- data.frame(stratum = NA, 
                                            death = NA,
                                            severe_disease = NA,
                                            infection = NA)
  for(i in 1:length(stratums)){
    I <- stratums[i]
    hold_efficacy <- filter(average_vacc_efficacy, stratum == I)
    hold_efficacy$death <- hold_efficacy$death *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$severe_disease <- hold_efficacy$severe_disease *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    hold_efficacy$infection <- hold_efficacy$infection *(sircovid:::sircovid_population("england")/sum(sircovid:::sircovid_population("england")))
    
    hold <- data.frame(stratum = I,
                       death = sum(hold_efficacy$death),
                       severe_disease = sum(hold_efficacy$severe_disease),
                       infection = sum(hold_efficacy$infection))
    
    average_vacc_efficacy_delta <- rbind(average_vacc_efficacy_delta, hold)
    
  }
  average_vacc_efficacy_delta <- average_vacc_efficacy_delta[-1,]
  
  
  #We now have all the VE, now we need to build the data frame we will plot
  #There will be six categories: 3*2 = (infection/hosp/death)*(fit/counterfactual)
  
  ###
  #Start by filtering out data to be as used in f2p1
  
  #First, we wish to only look at cumulative incidence
  keep_states <- c("first_dose", "second_dose", "waned_dose")
  
  f2p2_data <- filter(f2p2_data, state %in% keep_states)  
  
  #We will need to do some recalculating here too. As 1st dose cumulative should have
  #2nd dose removed from it, and waned dose should then be removed from second dose
  f2p2_data <- filter(f2p2_data, date > as.Date("2020-12-14"))
  f2p1_data_proc <- f2p2_data
  
  for(i in 1:length(f2p1_data_proc$state)){
    if(f2p1_data_proc$date[i] != as.Date("2020-12-15")){
      
      hold_state <- f2p1_data_proc$state[i]
      hold_date <- f2p1_data_proc$date[i]
      hold_type <- f2p1_data_proc$type[i]
      
      if(hold_state == "first_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "second_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }else if(hold_state == "second_dose"){
        hold <- filter(f2p1_data, date == hold_date & state == "waned_dose" & type == hold_type)
        f2p1_data_proc$value[i] <- f2p1_data_proc$value[i] - hold$value[1]
        
      }
      
    }
  }
  
  f2p1_data_proc$value <-  f2p1_data_proc$value/sum(sircovid:::sircovid_population("england"))
  #Now we make the new data frame that will be used for the fig, VE over time
  
  #JUST ALPHA HERE
  Alpha_VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "fit")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_alpha[,j])*(hold_population$value)))
      
      Alpha_VE_time_data_fit <- rbind(Alpha_VE_time_data_fit, hold)
    }
  }
  Alpha_VE_time_data_fit <- Alpha_VE_time_data_fit[-1,]
  
  
  #SAME BUT FOR COUNTERFACTUAL
  Alpha_VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "counterfactual")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_alpha[,j])*(hold_population$value)))
      
      Alpha_VE_time_data_counterfactual <- rbind(Alpha_VE_time_data_counterfactual, hold)
    }
  }
  Alpha_VE_time_data_counterfactual <- Alpha_VE_time_data_counterfactual[-1,]
  
  
  #Now for Delta
  Delta_VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "fit")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_delta[,j])*(hold_population$value)))
      
      Delta_VE_time_data_fit <- rbind(Delta_VE_time_data_fit, hold)
    }
  }
  Delta_VE_time_data_fit <- Delta_VE_time_data_fit[-1,]
  
  
  #SAME BUT FOR COUNTERFACTUAL
  Delta_VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA)
  outcomes <- c("death", "severe_disease", "infection")
  dates <- unique(f2p1_data_proc$date)
  for(i in 1:length(dates)){
    for(j in outcomes){
      I <- dates[i]
      #The population proportion in each stratum on date I:
      hold_population <- filter(f2p1_data_proc, date == I & type == "counterfactual")
      #3 values: 1dose, 2dose, waned
      hold <- data.frame(date = I,
                         outcome = j,
                         VE = sum((average_vacc_efficacy_delta[,j])*(hold_population$value)))
      
      Delta_VE_time_data_counterfactual <- rbind(Delta_VE_time_data_counterfactual, hold)
    }
  }
  Delta_VE_time_data_counterfactual <- Delta_VE_time_data_counterfactual[-1,]
  
  
  
  
  
  #Stick them together
  Delta_VE_time_data_fit$type <- "fit"
  Delta_VE_time_data_counterfactual$type <- "counterfactual"
  Alpha_VE_time_data_fit$type <- "fit"
  Alpha_VE_time_data_counterfactual$type <- "counterfactual"
  
  #Load in the Delta emergence date
  Delta_prop <- read.csv("Delta_emergence_data.csv")
  Delta_prop <- aggregate(Delta_prop$pos_mean, by = list(dates = Delta_prop$dates), FUN = mean)
  
  ####
  #And can directly extract the same for the counterfactual (in case it differs)
  
  n_strains <- filter(simulation_output$state, state %in% c("n_strain_1", "n_strain_2"))
  n_strains <- filter(n_strains, region == "england")
  n_strains <- filter(n_strains, quantile == "50%")
  n_strain_1 <- filter(n_strains, state == "n_strain_1")
  n_strain_2 <- filter(n_strains, state == "n_strain_2")
  
  n_strain_1$value <- c(0, diff(n_strain_1$value))
  n_strain_2$value <- c(0, diff(n_strain_2$value))
  n_strain_1 <- n_strain_1[-1,]
  n_strain_2 <- n_strain_2[-1,]
  
  n_strains <- rbind(n_strain_1, n_strain_2)
  
  simulation_delta_proportion <- data.frame(date = n_strain_1$date,
                                            prop_delta = n_strain_2$value/(n_strain_1$value + n_strain_2$value))
  
  simulation_delta_proportion <- filter(simulation_delta_proportion, date > as.Date("2021-03-07"))
  simulation_delta_proportion <- filter(simulation_delta_proportion, date < as.Date("2021-08-01"))
  ####
  
  #Build the final data frame
  VE_time_data_fit <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA, type = "fit")
  for(i in 1:length(dates)){
    I <- dates[i]
    if(I < as.Date("2021-03-08")){
      hold_Alpha <- filter(Alpha_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Alpha$VE[1]),
                                (hold_Alpha$VE[2]),
                                (hold_Alpha$VE[3])),
                         type = rep("fit", 3))
    } else if(I < as.Date("2021-08-01")){
      
      hold_Delta_prop <- Delta_prop$x[which(Delta_prop$dates == I)]/100
      hold_Alpha <- filter(Alpha_VE_time_data_fit, date == I )
      hold_Delta <- filter(Delta_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]*(hold_Delta_prop)) + (hold_Alpha$VE[1]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[2]*(hold_Delta_prop)) + (hold_Alpha$VE[2]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[3]*(hold_Delta_prop)) + (hold_Alpha$VE[3]*(1-hold_Delta_prop))),
                         type = rep("fit", 3))
      
    }else{
      
      hold_Delta <- filter(Delta_VE_time_data_fit, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]),
                                (hold_Delta$VE[2]),
                                (hold_Delta$VE[3])),
                         type = rep("fit", 3))
    }
    VE_time_data_fit <- rbind(VE_time_data_fit, hold)
    
  }
  VE_time_data_fit <- VE_time_data_fit[-1,]
  
  #REPEAT FOR COUNTERFACTUAL
  VE_time_data_counterfactual <- data.frame(date = as.Date("2000-12-12"), outcome = NA, VE = NA, type = "counterfactual")
  for(i in 1:length(dates)){
    I <- dates[i]
    if(I < as.Date("2021-03-08")){
      hold_Alpha <- filter(Alpha_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Alpha$VE[1]),
                                (hold_Alpha$VE[2]),
                                (hold_Alpha$VE[3])),
                         type = rep("counterfactual", 3))
    } else if(I < as.Date("2021-08-01")){
      
      hold_Delta_prop <- simulation_delta_proportion$prop_delta[which(simulation_delta_proportion$date == I)]
      hold_Alpha <- filter(Alpha_VE_time_data_counterfactual, date == I )
      hold_Delta <- filter(Delta_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]*(hold_Delta_prop)) + (hold_Alpha$VE[1]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[2]*(hold_Delta_prop)) + (hold_Alpha$VE[2]*(1-hold_Delta_prop)),
                                (hold_Delta$VE[3]*(hold_Delta_prop)) + (hold_Alpha$VE[3]*(1-hold_Delta_prop))),
                         type = rep("counterfactual", 3))
      
    }else{
      
      hold_Delta <- filter(Delta_VE_time_data_counterfactual, date == I )
      hold <- data.frame(date = rep(I,3),
                         outcome = c("death", "severe_disease", "infection"),
                         VE = c((hold_Delta$VE[1]),
                                (hold_Delta$VE[2]),
                                (hold_Delta$VE[3])),
                         type = rep("counterfactual", 3))
    }
    VE_time_data_counterfactual <- rbind(VE_time_data_counterfactual, hold)
    
  }
  VE_time_data_counterfactual <- VE_time_data_counterfactual[-1,]
  
  VE_time_data <- rbind(VE_time_data_fit, VE_time_data_counterfactual)
  
  #NOW we filter down for the stacked bars we want
  keep_dates <- seq(as.Date("2021/2/1"), by = "month", length.out = 9)
  VE_time_data_new <- filter(VE_time_data, date %in% keep_dates)
  
  #We can drop death
  VE_time_data_new <- filter(VE_time_data_new, outcome != "death")
 
  outcome.labs <- c("Infection", "Severe disease")
  names(outcome.labs) <- c("infection", "severe_disease")
  
  ggplot(data = VE_time_data_new, aes(x = date, y = VE)) +
    geom_col(data = filter(VE_time_data_new, type == "fit"), aes(x = date-6,group = type, fill = type), width = 10) +
    geom_col(data = filter(VE_time_data_new, type == "counterfactual"), aes(x = date+6, group = type, fill = type), width = 10) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,0.7)) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
    labs(fill = "Strategy", y = "Population vaccine protection", x = "Date") +
    geom_hline(yintercept = c(0.2, 0.4, 0.6), color = "black", alpha = 0.15) +
    facet_grid(rows = vars(outcome),
               labeller = labeller(outcome = outcome.labs)) +
    scale_fill_manual(values = c("#9BC362", "#B25A82"),
                      labels = c('3-week', "12-week")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="blank")
    ) 
  
   
}

#######
#EXTRA
#######

population_by_age_vacc <- function(simulation_output, combined){
  
  #Extract n_doses for each
  sim_n_doses <- simulation_output$n_doses
  fit_n_doses_raw <- combined$simulate$n_doses
  fit_n_doses_raw <- fit_n_doses_raw[,,8,]
  
  sim_n_doses <- filter(sim_n_doses, region == "england" & state != "booster_dose")
  sim_n_doses <- filter(sim_n_doses, date < as.Date("2021-09-14"))
  sim_n_doses <- subset(sim_n_doses, select=-c(vaccine_daily_doses,
                                               booster_daily_doses,
                                               beta_step,
                                               scenario, quantile,
                                               analysis, vaccine_status))
  
  n_doses_states <- dimnames(fit_n_doses_raw)[[2]]
  n_doses_ages <- unique(sim_n_doses$group)
  fit_n_doses <- data.frame(state = NA, region = NA, group = NA,
                            date = as.Date("2000-02-02") #placeholder!
                            , value = NA)
  #Age
  for(i in 1:19){
    #state
    for(j in 1:8){
      hold <- data.frame(state = rep(n_doses_states[j], 273), 
                         region = rep("england", 273), 
                         group = rep(n_doses_ages[i], 273),
                         date = sircovid::sircovid_date_as_date(combined$simulate$date), 
                         value = fit_n_doses_raw[i,j,])
      
      fit_n_doses <- rbind(fit_n_doses, hold)
    }
  }
  
  sim_n_doses$type <- rep("counterfactual", length(sim_n_doses$state))
  fit_n_doses$type <- rep("fit", length(fit_n_doses$state))
  #Remove the placeholder row
  fit_n_doses <- fit_n_doses[-1,]
  
  All_n_doses <- rbind(sim_n_doses, fit_n_doses)
  
  
  #Lastly, we aggregate into our specific age groups:
  #age_0 <- 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, (25% of CHW)
  #age_35 <- 35-39, 40-44, 45-49, 50-54, (50% of CHW)
  #age_55 <- 55-59, 60-64, 65-69, 70-74, (25% of CHW), (10% of CHR)
  #age_75 <- 75-79, 80+, (90% of CHR)
  
  All_n_doses$group <- as.factor(All_n_doses$group)
  init_levels <- levels(All_n_doses$group)
  new_levels <- c("age_0",  "age_0",  "age_0", "age_0", "age_0", "age_0", "age_0",
                  "age_35", "age_35", "age_35", "age_35",
                  "age_55", "age_55", "age_55", "age_55", "age_75", "age_75",
                  "CHW", "CHR")
  
  All_n_doses$group <- factor(new_levels[All_n_doses$group])

  #Aggregate by age now,
  
  All_n_doses <- aggregate(All_n_doses$value,  list(date = All_n_doses$date,
                                                 state = All_n_doses$state,
                                             group = All_n_doses$group,
                                                 type = All_n_doses$type),  sum)
  
  colnames(All_n_doses) <- c("date", "state", "group", "type", "value")
  
  All_n_doses <- filter(All_n_doses, state %in% c("first_dose", "second_dose", "waned_dose"))
  All_n_doses <- filter(All_n_doses, date > as.Date("2020-12-14"))
  
  #Now split up the CHW and CHR numbers
  dates_all <- unique(All_n_doses$date)
  states_all <- c("first_dose", "second_dose", "waned_dose")
  type_all <- unique(All_n_doses$type)
  
  for(i in 1:2){
    type_hold <- type_all[i]
    for(j in 1:length(dates_all)){
      dates_hold <- dates_all[j]
      for(k in 1:3){
        states_hold <- states_all[k]
        
        CHW_hold <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                      All_n_doses$date == dates_hold &
                                      All_n_doses$state == states_hold &
                                        All_n_doses$group == "CHW")]
        CHR_hold <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                              All_n_doses$date == dates_hold &
                                              All_n_doses$state == states_hold &
                                              All_n_doses$group == "CHR")]
        
        
        All_n_doses$value[which(All_n_doses$type == type_hold &
                                  All_n_doses$date == dates_hold &
                                  All_n_doses$state == states_hold &
                                  All_n_doses$group == "age_0")] <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                                                                              All_n_doses$date == dates_hold &
                                                                                              All_n_doses$state == states_hold &
                                                                                              All_n_doses$group == "age_0")] + 0.25*CHW_hold
        
        All_n_doses$value[which(All_n_doses$type == type_hold &
                                  All_n_doses$date == dates_hold &
                                  All_n_doses$state == states_hold &
                                  All_n_doses$group == "age_35")] <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                                                                              All_n_doses$date == dates_hold &
                                                                                              All_n_doses$state == states_hold &
                                                                                              All_n_doses$group == "age_35")] + 0.5*CHW_hold
        
        All_n_doses$value[which(All_n_doses$type == type_hold &
                                  All_n_doses$date == dates_hold &
                                  All_n_doses$state == states_hold &
                                  All_n_doses$group == "age_55")] <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                                                                              All_n_doses$date == dates_hold &
                                                                                              All_n_doses$state == states_hold &
                                                                                              All_n_doses$group == "age_55")] + 0.25*CHW_hold + 0.1*CHR_hold
        
        All_n_doses$value[which(All_n_doses$type == type_hold &
                                  All_n_doses$date == dates_hold &
                                  All_n_doses$state == states_hold &
                                  All_n_doses$group == "age_75")] <- All_n_doses$value[which(All_n_doses$type == type_hold &
                                                                                              All_n_doses$date == dates_hold &
                                                                                              All_n_doses$state == states_hold &
                                                                                              All_n_doses$group == "age_75")] + 0.9*CHR_hold
        
        
      }
    }
  }
  
  All_n_doses <- filter(All_n_doses, !(group %in% c("CHW", "CHR")))
  
  #Lastly, we need to remove second dose # from first dose #, and then waned dose # from second dose #
  
  All_n_doses$value[which(All_n_doses$state == "first_dose")] <- All_n_doses$value[which(All_n_doses$state == "first_dose")] - All_n_doses$value[which(All_n_doses$state == "second_dose")]
  All_n_doses$value[which(All_n_doses$state == "second_dose")] <- All_n_doses$value[which(All_n_doses$state == "second_dose")] - All_n_doses$value[which(All_n_doses$state == "waned_dose")]
  
  
  return(All_n_doses)
  
}









################################

prep_f2p3_data <- function(simulation_output, combined, p_C_and_H){
  
  #This plot will calculate population average risk of hospitalisation given infection.
  
  #p_C_and_H dimensions:
  #p_C: 1 x 19(age)
  #h_H:  7(region) x 1649(time) x 19(age) 
  #rel_p_hosp_if_sympt_VE: (19x4(strain)x4(vacc))
  #strain_rel_p_variant; (7(region) x 4 (strain))
  
  #We need to split the England population up by:
  #age
  #region
  #vacc status
  
  #Then can take average of these applied through the above parameters
  #factor in the Delta take-over
  #The only difference between the two scenarios is how the population is split
  #between vaccine strata
  
  #So first calculate, for each day, how many individuals (by age and region)
  #are in each vaccine strata. This is exactly what is done in f2p1
  
  
  
  
  #Extract n_doses for each
  sim_n_doses <- simulation_output$n_doses
  fit_n_doses_raw <- combined$simulate$n_doses
  
  sim_n_doses <- filter(sim_n_doses, region != "england" & state != "booster_dose")
  sim_n_doses <- filter(sim_n_doses, date < as.Date("2021-09-14"))
  sim_n_doses <- subset(sim_n_doses, select=-c(vaccine_daily_doses,
                                               booster_daily_doses,
                                               beta_step,
                                               scenario, quantile,
                                               analysis, vaccine_status))
  
  keep_states <- c("first_dose", "second_dose", "waned_dose")
  
  sim_n_doses <- filter(sim_n_doses, state %in% keep_states)  
  
  n_doses_states <- dimnames(fit_n_doses_raw)[[2]]
  n_doses_regions <- dimnames(fit_n_doses_raw)[[3]]
  n_doses_ages <- unique(sim_n_doses$group)
  fit_n_doses <- data.frame(state = NA, region = NA, group = NA,
                            date = as.Date("2000-02-02") #placeholder!
                            , value = NA)
  #Age
  for(i in 1:19){
    #Region
    for(k in 1:7){
    #state
    for(j in 1:3){
      hold <- data.frame(state = rep(n_doses_states[j], 273), 
                         region = rep(n_doses_regions[k], 273), 
                         group = rep(n_doses_ages[i], 273),
                         date = sircovid::sircovid_date_as_date(combined$simulate$date), 
                         value = fit_n_doses_raw[i,j,k,])
      
      fit_n_doses <- rbind(fit_n_doses, hold)
    }
    }
  }
  
  
  sim_n_doses$type <- rep("counterfactual", length(sim_n_doses$state))
  fit_n_doses$type <- rep("fit", length(fit_n_doses$state))
  #Remove the placeholder row
  fit_n_doses <- fit_n_doses[-1,]
  
  All_n_doses <- rbind(sim_n_doses, fit_n_doses)
  
  
  #we only plan to use the first of each month
  keep_dates <- seq(as.Date("2021/2/1"), by = "month", length.out = 8)
  
  All_n_doses <- filter(All_n_doses, date %in% keep_dates)
  
  #We now need to remove second dose # from first dose # and then waned from second dose
  
  for(i in 1:length(keep_dates)){
    for(j in 1:length(n_doses_ages)){
      for(k in n_doses_regions[1:7]){
        for(l in c("fit", "counterfactual")){
          
          All_n_doses[which(All_n_doses$date == keep_dates[i] &
                            All_n_doses$group == n_doses_ages[j] &
                            All_n_doses$region == k &
                            All_n_doses$type == l &
                            All_n_doses$state == "first_dose"),]$value <- All_n_doses[which(All_n_doses$date == keep_dates[i] &
                                                                                              All_n_doses$group == n_doses_ages[j] &
                                                                                              All_n_doses$region == k &
                                                                                              All_n_doses$type == l &
                                                                                              All_n_doses$state == "first_dose"),]$value - All_n_doses[which(All_n_doses$date == keep_dates[i] &
                                                                                                                                                               All_n_doses$group == n_doses_ages[j] &
                                                                                                                                                               All_n_doses$region == k &
                                                                                                                                                               All_n_doses$type == l &
                                                                                                                                                               All_n_doses$state == "second_dose"),]$value
          All_n_doses[which(All_n_doses$date == keep_dates[i] &
                              All_n_doses$group == n_doses_ages[j] &
                              All_n_doses$region == k &
                              All_n_doses$type == l &
                              All_n_doses$state == "second_dose"),]$value <- All_n_doses[which(All_n_doses$date == keep_dates[i] &
                                                                                                 All_n_doses$group == n_doses_ages[j] &
                                                                                                 All_n_doses$region == k &
                                                                                                 All_n_doses$type == l &
                                                                                                 All_n_doses$state == "second_dose"),]$value - All_n_doses[which(All_n_doses$date == keep_dates[i] &
                                                                                                                                                                  All_n_doses$group == n_doses_ages[j] &
                                                                                                                                                                  All_n_doses$region == k &
                                                                                                                                                                  All_n_doses$type == l &
                                                                                                                                                                  All_n_doses$state == "waned_dose"),]$value
          
        }
      }
    }
  }
  
  
  #We need to now add in the unvaccinated population.
  #For each day / age / region / type
  unvacc_data <- data.frame(state = "unvaccinated",
                            region = "",
                            group = "",
                            date = keep_dates[1],
                            value = 0,
                            type = "")
  for(i in 1:length(keep_dates)){
    for(j in 1:length(n_doses_ages)){
      for(k in n_doses_regions[1:7]){
        for(l in c("fit", "counterfactual")){
          
          All_n_doses %>%
            filter(date == keep_dates[i]) %>%
            filter(group == n_doses_ages[j]) %>%
            filter(region == k) %>%
            filter(type == l) -> hold
          
          #Add in CHR CHW to population tally
          population <- sircovid::lancelot_parameters(1, k)$N_tot
          
          
          #first dose number
          first_dose_hold <- filter(hold, state == "first_dose")
          second_dose_hold <- filter(hold, state == "second_dose")
          waned_dose_hold <- filter(hold, state == "waned_dose")
          
          unvacc_hold <- data.frame(state = "unvaccinated",
                                    region = k,
                                    group = n_doses_ages[j],
                                    date = keep_dates[i],
                                    value = population[j] - first_dose_hold$value - second_dose_hold$value - waned_dose_hold$value,
                                    type = l)
          
          unvacc_data <- rbind(unvacc_data, unvacc_hold)
          
        }
      }
    }
  }
  
  unvacc_data <- unvacc_data[-1,]
  
  All_n_doses <- rbind(All_n_doses, unvacc_data)
  
  
  #Load in the Delta emergence date
  Delta_prop <- read.csv("Delta_emergence_data.csv")
  Delta_prop <- aggregate(Delta_prop$pos_mean, by = list(dates = Delta_prop$dates), FUN = mean)
  Delta_prop$x <- Delta_prop$x/100
  colnames(Delta_prop) <- c("date", "prop_delta")
  Delta_prop$type <- "fit"
  Delta_prop$date <- as.Date(Delta_prop$date)
  Delta_prop <- filter(Delta_prop, date %in% keep_dates)
  ####
  #And can directly extract the same for the counterfactual (in case it differs)
  
  n_strains <- filter(simulation_output$state, state %in% c("n_strain_1", "n_strain_2"))
  n_strains <- filter(n_strains, region == "england")
  n_strains <- filter(n_strains, quantile == "50%")
  n_strain_1 <- filter(n_strains, state == "n_strain_1")
  n_strain_2 <- filter(n_strains, state == "n_strain_2")
  
  n_strain_1$value <- c(0, diff(n_strain_1$value))
  n_strain_2$value <- c(0, diff(n_strain_2$value))
  n_strain_1 <- n_strain_1[-1,]
  n_strain_2 <- n_strain_2[-1,]
  
  n_strains <- rbind(n_strain_1, n_strain_2)
  
  simulation_delta_proportion <- data.frame(date = n_strain_1$date,
                                            prop_delta = n_strain_2$value/(n_strain_1$value + n_strain_2$value))
  
  simulation_delta_proportion <- filter(simulation_delta_proportion, date > as.Date("2021-03-07"))
  simulation_delta_proportion <- filter(simulation_delta_proportion, date < as.Date("2021-08-01"))
  
  simulation_delta_proportion$type <- "counterfactual"
  simulation_delta_proportion <- filter(simulation_delta_proportion, date %in% keep_dates)
  
  sim_add <- data.frame(date = keep_dates[c(1,2,7,8)],
                        prop_delta = c(0,0,1,1),
                        type = rep("counterfactual"))
  simulation_delta_proportion <- rbind(simulation_delta_proportion, sim_add)
  
  fit_add <- data.frame(date = keep_dates[c(1,2,7,8)],
                        prop_delta = c(0,0,1,1),
                        type = rep("fit"))
  Delta_prop <- rbind(Delta_prop, fit_add)
  Delta_prop <- rbind(Delta_prop, simulation_delta_proportion)
  ####
  
  vacc_states <- c("unvaccinated", "first_dose", "second_dose", "waned_dose")
  
  #We now add a new column that contains the p_C * p_H value for each row
  #p_C_and_H dimensions:
  #p_C: 1 x 19(age)
  #h_H:  7(region) x 1649(time) x 19(age) 
  #rel_p_hosp_if_sympt_VE: (19x4(strain)x4(vacc))
  #strain_rel_p_variant; (7(region) x 4 (strain))
  p_C_H <- rep(0, length(All_n_doses$state))
  each_population <- rep(0, length(All_n_doses$state))
  
  for(i in 1:length(All_n_doses$state)){
    age_hold <- All_n_doses$group[i]
    region_hold <- All_n_doses$region[i]
    vacc_hold <- All_n_doses$state[i]
    time_hold <- All_n_doses$date[i]
    type_hold <- All_n_doses$type[i]
    sircovid_time_hold <- sircovid::sircovid_date(time_hold)*4
    if(sircovid_time_hold > 1649){
      sircovid_time_hold <- 1649
    }
    pop_hold <- sircovid::lancelot_parameters(1, region_hold)$N_tot
    each_population[i] <- pop_hold[which(n_doses_ages == age_hold)]
    
    p_C <- p_C_and_H$p_C[which(n_doses_ages == age_hold)]
    p_H <- p_C_and_H$h_H[which(n_doses_regions == region_hold),
                                         sircovid_time_hold,
                                         which(n_doses_ages == age_hold)]
    
    #Calculate the Alpha and the Delta strain values
    #Sadly we don't have reinfections listed here, so treat each one as a first infection
    
    rel_p_alpha <- p_C_and_H$rel_p_hosp_if_sympt_VE[which(n_doses_ages == age_hold),
                                                    1,
                                                    which(vacc_states == vacc_hold)]
    
    rel_p_delta <- p_C_and_H$rel_p_hosp_if_sympt_VE[which(n_doses_ages == age_hold),
                                                    2,
                                                    which(vacc_states == vacc_hold)]
    
    prop_data <- filter(Delta_prop, type == type_hold)
    prop_data <- filter(prop_data, date == time_hold)
    
    rel_p <- rel_p_delta*prop_data$prop_delta[1] + rel_p_alpha*(1-prop_data$prop_delta[1])
    
    p_H <- p_H*rel_p
    
    strain_rel_alpha <- p_C_and_H$strain_rel_p_variant[which(n_doses_regions == region_hold),1]
    strain_rel_delta <- p_C_and_H$strain_rel_p_variant[which(n_doses_regions == region_hold),2]
    
    strain_rel <- strain_rel_delta*prop_data$prop_delta[1] + strain_rel_alpha*(1-prop_data$prop_delta[1])
    
    p_H <- p_H*strain_rel
    p_H <- min(1, p_H)
    
    p_C_H[i] <- p_C*p_H
  }
  
  All_n_doses$p_C_and_H <- p_C_H
  All_n_doses$population <- each_population
  
  #Lastly, we generate the actual data frame we want. 
  #16 rows, 8 for each scenario, 8 dates. We need to get the mean p_C_H
  #over age and region
  f2p3_data <- data.frame(date = c(keep_dates, keep_dates),
                          type = c(rep("fit", 8), rep("counterfactual", 8)),
                          p_C_H = rep(0,16))
  
  for(i in 1:16){
    date_hold <- f2p3_data$date[i]
    type_hold <- f2p3_data$type[i]
    
    doses_data <- filter(All_n_doses, (date == date_hold & type == type_hold ))
    
    total_risk <- sum(doses_data$p_C_and_H*doses_data$value)/sum(doses_data$value)
  f2p3_data$p_C_H[i] <- total_risk
      
  }
  
  return(f2p3_data)
  
}






############
plot_f2p3 <- function(f2p3_data){
  
  ggplot(data = f2p3_data,
         aes(x = date, y = p_C_H)) +
    geom_col(data = filter(f2p3_data, type == "fit"), aes(x = date-6,group = type, fill = type), width = 10) +
    geom_col(data = filter(f2p3_data, type == "counterfactual"), aes(x = date+6, group = type, fill = type), width = 10) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,0.09),
                       breaks=c(0.02, 0.04, 0.06, 0.08)) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
    labs(fill = "Strategy", y = "Population average risk of \nhospitalisation given infection", x = "Date") +
    geom_hline(yintercept = c(0.02, 0.04, 0.06, 0.08), color = "black", alpha = 0.15) +
    #facet_grid(rows = vars(outcome),
    #           labeller = labeller(outcome = outcome.labs)) +
    scale_fill_manual(breaks = c("fit", "counterfactual"),
                      values = c("#B25A82", "#9BC362"),
                      labels = c("12-week", '3-week')) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(
            color="black", linetype="blank")
    ) 
  
}










