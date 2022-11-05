

plot_Alpha <- function(f4_data_alpha){
  
  #First, provide the data being fit to:
  
  tau <- 24
  
  ###############################################
  
  
  y_Mild_AZ_ALPHA <- c(0.718, 0.824, 0.762, NA, NA)
  y_Mild_AZ_ALPHA_lower <- c(0.662, 0.796, 0.498, NA, NA)
  y_Mild_AZ_ALPHA_upper <- c(0.765, 0.847, 0.887, NA, NA)
  y_Mild_AZ_ALPHA_CI_range <- y_Mild_AZ_ALPHA_upper - y_Mild_AZ_ALPHA_lower
  
  y_Mild_PF_ALPHA <- c(0.902, 0.949, 0.948, NA, NA)
  y_Mild_PF_ALPHA_lower <- c(0.869, 0.936, 0.884, NA, NA)
  y_Mild_PF_ALPHA_upper <- c(0.927, 0.959, 0.977, NA, NA)
  y_Mild_PF_ALPHA_CI_range <- y_Mild_PF_ALPHA_upper - y_Mild_PF_ALPHA_lower
  
  y_Severe_AZ_ALPHA <- c(0.896, 0.951, NA, NA, NA)
  y_Severe_AZ_ALPHA_lower <- c(0.674, 0.867, NA, NA, NA)
  y_Severe_AZ_ALPHA_upper <- c(0.967, 0.982, NA, NA, NA)
  y_Severe_AZ_ALPHA_CI_range <- y_Severe_AZ_ALPHA_upper - y_Severe_AZ_ALPHA_lower
  
  y_Severe_PF_ALPHA <- c(0.964, 0.977, NA, NA, NA)
  y_Severe_PF_ALPHA_lower <- c(0.742, 0.908, NA, NA, NA)
  y_Severe_PF_ALPHA_upper <- c(0.995, 0.994, NA, NA, NA)
  y_Severe_PF_ALPHA_CI_range <- y_Severe_PF_ALPHA_upper - y_Severe_PF_ALPHA_lower
  
  y_Death_AZ_ALPHA <- c(NA, 0.961, NA, NA, NA, NA)
  y_Death_AZ_ALPHA_lower <- c(NA, 0.89, NA, NA, NA, NA)
  y_Death_AZ_ALPHA_upper <- c(NA, 0.986, NA, NA, NA, NA)
  y_Death_AZ_ALPHA_CI_range <- y_Death_AZ_ALPHA_upper - y_Death_AZ_ALPHA_lower
  
  y_Death_PF_ALPHA <- c(NA, 0.961, NA, NA, NA, NA)
  y_Death_PF_ALPHA_lower <- c(NA, 0.89, NA, NA, NA, NA)
  y_Death_PF_ALPHA_upper <- c(NA, 0.986, NA, NA, NA, NA)
  y_Death_PF_ALPHA_CI_range <- y_Death_PF_ALPHA_upper - y_Death_PF_ALPHA_lower
  
  Andrews_Data <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                             timeframe = c("2-9",
                                           "10-14",
                                           "15-19",
                                           "20+"),
                             sensitivity = rep("data", 4),
                             outcome = rep("mild", 4),
                             vaccine = rep("AZ",4),
                             variant = rep("Alpha", 4),
                             median = y_Mild_AZ_ALPHA[2:5],
                             lower = y_Mild_AZ_ALPHA_lower[2:5],
                             upper = y_Mild_AZ_ALPHA_upper[2:5])
  
  Andrews_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("mild", 4),
                                  vaccine = rep("PF",4),
                                  variant = rep("Alpha", 4),
                                  median = y_Mild_PF_ALPHA[2:5],
                                  lower = y_Mild_PF_ALPHA_lower[2:5],
                                  upper = y_Mild_PF_ALPHA_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Severe AZ
  Andrews_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("severe", 4),
                                  vaccine = rep("AZ",4),
                                  variant = rep("Alpha", 4),
                                  median = y_Severe_AZ_ALPHA[2:5],
                                  lower = y_Severe_AZ_ALPHA_lower[2:5],
                                  upper = y_Severe_AZ_ALPHA_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Severe PF
  Andrews_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("severe", 4),
                                  vaccine = rep("PF",4),
                                  variant = rep("Alpha", 4),
                                  median = y_Severe_PF_ALPHA[2:5],
                                  lower = y_Severe_PF_ALPHA_lower[2:5],
                                  upper = y_Severe_PF_ALPHA_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Death AZ
  Andrews_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("death", 4),
                                  vaccine = rep("AZ",4),
                                  variant = rep("Alpha", 4),
                                  median = y_Death_AZ_ALPHA[2:5],
                                  lower = y_Death_AZ_ALPHA_lower[2:5],
                                  upper = y_Death_AZ_ALPHA_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Death PF
  Andrews_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("death", 4),
                                  vaccine = rep("PF",4),
                                  variant = rep("Alpha", 4),
                                  median = y_Death_PF_ALPHA[2:5],
                                  lower = y_Death_PF_ALPHA_lower[2:5],
                                  upper = y_Death_PF_ALPHA_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  rm(Andrews_Data_hold)
  
  #Now build a data frame with our simulation equivalents;
  ###########################################################################
  ###########################################################################
  
  #Mild AZ
  ###
  #baseline
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)])))
  
  Sim_Data <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("baseline", 4),
                         outcome = rep("mild", 4),
                         vaccine = rep("AZ",4),
                         variant = rep("Alpha", 4),
                         median = c(mean(waned_VE[(14-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(202-13)])),
                         lower = rep(NA, 4),
                         upper = rep(NA, 4))
  
  
  
  ################
  #PF MILD DELTA
  ###
  #baseline
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                              timeframe = c("2-9", #Days 14-69
                                            "10-14", #Days 70-104
                                            "15-19", #Days 105-139
                                            "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("mild", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Alpha", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  #Severe AZ
  ###
  #baseline
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Alpha", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  
  
  
  #Severe PF
  ###
  #baseline
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Alpha", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  
  
  
  #Death PF
  ###
  #baseline
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Alpha", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  
  
  
  #Death AZ
  ###
  #baseline
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Alpha", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  ############################################
  ## PLOTTING
  ############################################
  
  
  #AZ Death
  
  filter(Andrews_Data, label == "Alpha_Death_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Death_AZ")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe == "2-9")
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4), 
               stroke = 1.5,
               #shape = 18, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - AZ - Death") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792","#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE)) +
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> AZ_death
  
  #PF Death
  
  filter(Andrews_Data, label == "Alpha_Death_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Death_PF")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe == "2-9")
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4), 
               stroke = 1.5, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    #ylim(c(0.6,1)) +
    ggtitle("Alpha - PF - Death") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792", "#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> PF_death
  
  
  #AZ Severe
  
  filter(Andrews_Data, label == "Alpha_Severe_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Severe_AZ")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe == "2-9")
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4), 
               stroke = 1.5, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    #ylim(c(0.6,1)) +
    ggtitle("Alpha - AZ - Severe") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792", "#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> AZ_severe
  
  #PF Severe
  
  filter(Andrews_Data, label == "Alpha_Severe_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Severe_PF")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe == "2-9")
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4),
               stroke = 1.5, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    #ylim(c(0.6,1)) +
    ggtitle("Alpha - PF - Severe") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792", "#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> PF_severe
  
  
  #AZ Mild
  
  filter(Andrews_Data, label == "Alpha_Mild_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Mild_AZ")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe %in% c("2-9", "10-14"))
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4), 
               stroke = 1.5, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    #ylim(c(0.2,1)) +
    ggtitle("Alpha - AZ - Mild") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792", "#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = c(0.2, 0.6, 1),
                       limits = c(0.2, 1)) +
    guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> AZ_mild
  
  
  #PF Mild
  
  filter(Andrews_Data, label == "Alpha_Mild_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline")
  d_simulated <- filter(Sim_Data, label == "Alpha_Mild_PF")
  d_simulated <- filter(d_simulated, sensitivity == "baseline")
  d_simulated <- filter(d_simulated, timeframe %in% c("2-9", "10-14"))
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(data_hold, aes(group = sensitivity)) +
    geom_errorbar(aes(x= timeframe, ymin = lower, ymax = upper,
                      color = sensitivity, group = sensitivity),
                  #position = position_dodge(width=0.4),
                  #color = "#841744",
                  width = 0.3,
                  alpha = 1) +
    geom_point(aes(x = timeframe, y = median, color = sensitivity,
                   shape = sensitivity, group = sensitivity),
               #position = position_dodge(width=0.4), 
               stroke = 1.5, 
               size = 4,
               alpha = 1) +
    theme_classic() +
    
    xlab("Weeks since second dose") + 
    #ylim(c(0.2,1)) +
    ggtitle("Alpha - PF - Mild") +
    ylab("VE") +
    labs(colour="") +
    scale_color_manual(values=c(#"#7be8bb", 
      "#1a9792", "#841744"),
      labels = c("Model", "Data")) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(1, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Model", "Data"),
                       values= c(5, 16)) +
    scale_y_continuous(breaks = c(0.2, 0.6, 1),
                       limits = c(0.2, 1)) +
    guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) -> PF_mild
  
  prow <- plot_grid(
    AZ_death + theme(legend.position="none") + xlab(NULL),
    PF_death + theme(legend.position="none") + xlab(NULL) + ylab(NULL),
    AZ_severe + theme(legend.position="none") + xlab(NULL),
    PF_severe + theme(legend.position="none") + xlab(NULL) + ylab(NULL),
    AZ_mild + theme(legend.position="none"),
    PF_mild + theme(legend.position="none") + ylab(NULL),
    align = 'vh',
    labels = c("A", "B", "C", "D", "E", "F"),
    hjust = -1,
    nrow = 3
  )
  
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    AZ_death + theme(legend.box.margin = margin(0, 0, 0, 4))
  )
  
  # add the legend to the row we made earlier. 
  plot_grid(prow, legend, rel_widths = c(3, 0.8))
  
  
}



#These functions take in a mean week waning, and VE parameters, they output the disease trajectory
#for 230 days

death_plotting <- function(par){
  
  tau <- par[1]
  two_dose_VE_death <- par[2]
  two_dose_VE_severe <- par[3]
  two_dose_VE_mild <- par[4]
  waned_VE_death <- par[5]
  waned_VE_severe <- par[6]
  waned_VE_mild <- par[7]
  
  
  days <- 1:230
  proportion_waned <- pexp(days, rate = 1/(tau*7))
  population_death_VE <- (waned_VE_death*proportion_waned) + (two_dose_VE_death*(1-proportion_waned))
  
  
  return(population_death_VE)
  
}

severe_plotting <- function(par){
  
  tau <- par[1]
  two_dose_VE_death <- par[2]
  two_dose_VE_severe <- par[3]
  two_dose_VE_mild <- par[4]
  waned_VE_death <- par[5]
  waned_VE_severe <- par[6]
  waned_VE_mild <- par[7]
  
  
  days <- 1:230
  proportion_waned <- pexp(days, rate = 1/(tau*7))
  population_severe_VE <- (waned_VE_severe*proportion_waned) + (two_dose_VE_severe*(1-proportion_waned))
  
  return(population_severe_VE)
  
}

mild_plotting <- function(par){
  
  tau <- par[1]
  two_dose_VE_death <- par[2]
  two_dose_VE_severe <- par[3]
  two_dose_VE_mild <- par[4]
  waned_VE_death <- par[5]
  waned_VE_severe <- par[6]
  waned_VE_mild <- par[7]
  
  
  days <- 1:230
  proportion_waned <- pexp(days, rate = 1/(tau*7))
  population_mild_VE <- (waned_VE_mild*proportion_waned) + (two_dose_VE_mild*(1-proportion_waned))
  
  return(population_mild_VE)
  
}


