

plot_f6a <- function(f4_data_delta){
  
  #First, provide the data being fit to:
  
  tau <- 24
  
  ###############################################
  # Delta symptomatic taken from 10.1056/NEJMoa2119451
  #All others taken from 10.1056/NEJMoa2115481
  #The six VE estimates correspond to time periods:
  #2-4 weeks (treated as Days 14-34)
  #5-9 weeks (treated as Days 35-69)
  #10- 14 weeks (treated as Days 70-104)
  #15-19 weeks (treated as Days 105-139)
  #20-24 weeks (treated as Days 140-174)
  #25+ weeks  (treated as Days 175-209)
  y_Mild_AZ <- c(0.828, 0.765, 0.692, 0.536, 0.474, 0.435)
  y_Mild_AZ_lower <- c(0.745, 0.703, 0.647, 0.516, 0.462, 0.424)
  y_Mild_AZ_upper <- c(0.884, 0.815, 0.731, 0.555, 0.485, 0.445)
  y_Mild_AZ_CI_range <- y_Mild_AZ_upper - y_Mild_AZ_lower
  
  y_Mild_PF <- c(0.909, 0.855, 0.787, 0.744, 0.674, 0.627)
  y_Mild_PF_lower <- c(0.896, 0.845, 0.78, 0.738, 0.665, 0.616)
  y_Mild_PF_upper <- c(0.92, 0.865, 0.794, 0.749, 0.682, 0.637)
  y_Mild_PF_CI_range <- y_Mild_PF_upper - y_Mild_PF_lower
  
  #The five VE estimates correspond to time periods:
  #<2 weeks (treated as Days 1-13)
  #2-9 weeks (treated as Days 14-69)
  #10 - 14 weeks (treated as Days 70-104)
  #15-19 weeks (treated as Days 105-139)
  #20+ weeks (treated as Days 140-202)
  y_Severe_AZ <- c(0.94, 0.952, 0.921, 0.874, 0.800)
  y_Severe_AZ_lower <- c(0.913, 0.947, 0.913, 0.861, 0.768)
  y_Severe_AZ_upper <- c(0.958, 0.957, 0.927, 0.886, 0.827)
  y_Severe_AZ_CI_range <- y_Severe_AZ_upper - y_Severe_AZ_lower
  
  y_Severe_PF <- c(0.994, 0.987, 0.968, 0.949, 0.917)
  y_Severe_PF_lower <- c(0.977, 0.983, 0.963, 0.941, 0.902)
  y_Severe_PF_upper <- c(0.999, 0.99, 0.973, 0.955, 0.93)
  y_Severe_PF_CI_range <- y_Severe_PF_upper - y_Severe_PF_lower
  
  y_Death_AZ <- c(NA, 0.95, 0.937, 0.901, 0.848)
  y_Death_AZ_lower <- c(NA, 0.931, 0.918, 0.869, 0.762)
  y_Death_AZ_upper <- c(NA, 0.964, 0.952, 0.926, 0.903)
  y_Death_AZ_CI_range <- y_Death_AZ_upper - y_Death_AZ_lower
  
  y_Death_PF <- c(NA, 0.985, 0.96, 0.945, 0.919)
  y_Death_PF_lower <- c(NA, 0.965, 0.942, 0.925, 0.885)
  y_Death_PF_upper <- c(NA, 0.993, 0.972, 0.96, 0.943)
  y_Death_PF_CI_range <- y_Death_PF_upper - y_Death_PF_lower
  
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
  
  Andrews_Data <- data.frame(label = rep("Delta_Mild_AZ", 6),
                             timeframe = c("2-4",
                                           "5-9",
                                           "10-14",
                                           "15-19",
                                           "20-24",
                                           "25+"),
                             sensitivity = rep("data", 6),
                             outcome = rep("mild", 6),
                             vaccine = rep("AZ",6),
                             variant = rep("Delta", 6),
                             median = y_Mild_AZ,
                             lower = y_Mild_AZ_lower,
                             upper = y_Mild_AZ_upper)
  
  Andrews_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                                  timeframe = c("2-4",
                                                "5-9",
                                                "10-14",
                                                "15-19",
                                                "20-24",
                                                "25+"),
                                  sensitivity = rep("data", 6),
                                  outcome = rep("mild", 6),
                                  vaccine = rep("PF",6),
                                  variant = rep("Delta", 6),
                                  median = y_Mild_PF,
                                  lower = y_Mild_PF_lower,
                                  upper = y_Mild_PF_upper)
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Severe AZ
  Andrews_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("severe", 4),
                                  vaccine = rep("AZ",4),
                                  variant = rep("Delta", 4),
                                  median = y_Severe_AZ[2:5],
                                  lower = y_Severe_AZ_lower[2:5],
                                  upper = y_Severe_AZ_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Severe PF
  Andrews_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("severe", 4),
                                  vaccine = rep("PF",4),
                                  variant = rep("Delta", 4),
                                  median = y_Severe_PF[2:5],
                                  lower = y_Severe_PF_lower[2:5],
                                  upper = y_Severe_PF_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Death AZ
  Andrews_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("death", 4),
                                  vaccine = rep("AZ",4),
                                  variant = rep("Delta", 4),
                                  median = y_Death_AZ[2:5],
                                  lower = y_Death_AZ_lower[2:5],
                                  upper = y_Death_AZ_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  ####################
  #Death PF
  Andrews_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                                  timeframe = c("2-9",
                                                "10-14",
                                                "15-19",
                                                "20+"),
                                  sensitivity = rep("data", 4),
                                  outcome = rep("death", 4),
                                  vaccine = rep("PF",4),
                                  variant = rep("Delta", 4),
                                  median = y_Death_PF[2:5],
                                  lower = y_Death_PF_lower[2:5],
                                  upper = y_Death_PF_upper[2:5])
  
  Andrews_Data <- rbind(Andrews_Data, Andrews_Data_hold)
  
  rm(Andrews_Data_hold)
  
  #Now build a data frame with our simulation equivalents;
  ###########################################################################
  ###########################################################################
  
  #Mild AZ
  ###
  #baseline
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data <- data.frame(label = rep("Delta_Mild_AZ", 6),
                         timeframe = c("2-4",#(treated as Days 14-34)
                                       "5-9", #(treated as Days 35-69)
                                       "10-14", #(treated as Days 70-104)
                                       "15-19", #(treated as Days 105-139)
                                       "20-24", #(treated as Days 140-174)
                                       "25+"), #(treated as Days 175-209)
                         sensitivity = rep("baseline", 6),
                         outcome = rep("mild", 6),
                         vaccine = rep("AZ",6),
                         variant = rep("Delta", 6),
                         median = c(mean(waned_VE[(14-13):(34-13)]),
                                    mean(waned_VE[(35-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(174-13)]),
                                    mean(waned_VE[(175-13):(209-13)])),
                         lower = rep(NA, 6),
                         upper = rep(NA, 6))
  
  
  #sensitivity #4
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#4 - all upper", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#5 - all lower", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#6 - waned upper", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#7 - waned lower", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#2 - reduced mild", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#3 - 0.1 reduction all", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- mild_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_AZ", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#8 - 33 week waning", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("AZ",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  ################
  #PF MILD DELTA
  ###
  #baseline
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("baseline", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #4
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#4 - all upper", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#5 - all lower", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#6 - waned upper", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#7 - waned lower", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#2 - reduced mild", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#3 - 0.1 reduction all", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- mild_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Mild_PF", 6),
                              timeframe = c("2-4",#(treated as Days 14-34)
                                            "5-9", #(treated as Days 35-69)
                                            "10-14", #(treated as Days 70-104)
                                            "15-19", #(treated as Days 105-139)
                                            "20-24", #(treated as Days 140-174)
                                            "25+"), #(treated as Days 175-209)
                              sensitivity = rep("#8 - 33 week waning", 6),
                              outcome = rep("mild", 6),
                              vaccine = rep("PF",6),
                              variant = rep("Delta", 6),
                              median = c(mean(waned_VE[(14-13):(34-13)]),
                                         mean(waned_VE[(35-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(174-13)]),
                                         mean(waned_VE[(175-13):(209-13)])),
                              lower = rep(NA, 6),
                              upper = rep(NA, 6))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  
  
  
  #Severe AZ
  ###
  #baseline
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #4
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#4 - all upper", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#5 - all lower", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#6 - waned upper", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#7 - waned lower", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#2 - reduced mild", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#3 - 0.1 reduction all", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- severe_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#8 - 33 week waning", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
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
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #4
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#4 - all upper", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#5 - all lower", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#6 - waned upper", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#7 - waned lower", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#2 - reduced mild", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#3 - 0.1 reduction all", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- severe_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Severe_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#8 - 33 week waning", 4),
                              outcome = rep("severe", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
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
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #4
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#4 - all upper", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#5 - all lower", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#6 - waned upper", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#7 - waned lower", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#2 - reduced mild", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#3 - 0.1 reduction all", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- death_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_PF", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#8 - 33 week waning", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("PF",4),
                              variant = rep("Delta", 4),
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
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ), 
                              sensitivity = rep("baseline", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #4
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#4 - all upper", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#5 - all lower", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #6 - waned upper
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_6_higher_wane[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#6 - waned upper", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #7
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_7_lower_waned[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#7 - waned lower", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#2 - reduced mild", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_delta$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#3 - 0.1 reduction all", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #8
  waned_VE <- death_plotting(c(33, as.numeric(f4_data_delta$central[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Delta_Death_AZ", 4),
                              timeframe = c(
                                "2-9", #Days 14-69
                                "10-14", #Days 70-104
                                "15-19", #Days 105-139
                                "20+"#Days 140-202
                              ),
                              sensitivity = rep("#8 - 33 week waning", 4),
                              outcome = rep("death", 4),
                              vaccine = rep("AZ",4),
                              variant = rep("Delta", 4),
                              median = c(mean(waned_VE[(14-13):(69-13)]),
                                         mean(waned_VE[(70-13):(104-13)]),
                                         mean(waned_VE[(105-13):(139-13)]),
                                         mean(waned_VE[(140-13):(202-13)])),
                              lower = rep(NA, 4),
                              upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  Sim_Data$sensitivity <- factor(Sim_Data$sensitivity,
                                    levels = c("baseline", "#2 - reduced mild",
                                               "#3 - 0.1 reduction all", "#4 - all upper",        
                                               "#5 - all lower",   "#6 - waned upper",      
                                               "#7 - waned lower", "#8 - 33 week waning" ))
  ############################################
  ## PLOTTING
  ############################################
  
  
  #AZ Death
  
  filter(Andrews_Data, label == "Delta_Death_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower", "#6 - waned upper",
                          "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Death_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
               shape = sensitivity, group = sensitivity, shape = sensitivity),
              stroke = 1.5,
               color = "#841744",
               size = 2,
              key_glyph = "pointrange" ) +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3
                  ) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - AZ - Death") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    # scale_color_lancet(labels = c("#1 - Baseline",
    #                               "#3 - Absolute decrease",
    #                               "#4 - Relative increase",
    #                               "#5 - Relative decrease",
    #                               "#6 - Waned increased",
    #                               "#7 - Waned decreased"
    #                               )) +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
                                  ),
                       labels = c("#1 - Baseline",
                                  "#3 - Absolute decrease",
                                  "#4 - Relative increase",
                                  "#5 - Relative decrease",
                                  "#6 - Waned increased",
                                  "#7 - Waned decreased"
                       )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_death
  
  #PF Death
  
  filter(Andrews_Data, label == "Delta_Death_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower", "#6 - waned upper",
                          "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Death_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - PF - Death") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease",
               "#6 - Waned increased",
               "#7 - Waned decreased"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_death
  
  #AZ Severe
  
  filter(Andrews_Data, label == "Delta_Severe_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower", "#6 - waned upper",
                          "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Severe_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - AZ - Severe Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease",
               "#6 - Waned increased",
               "#7 - Waned decreased"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_severe
  
  #PF Severe
  
  filter(Andrews_Data, label == "Delta_Severe_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- sensitivities_keep <- c("baseline", "#4 - all upper",
                                                "#3 - 0.1 reduction all",
                                                "#5 - all lower", "#6 - waned upper",
                                                "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Severe_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - PF - Severe Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease",
               "#6 - Waned increased",
               "#7 - Waned decreased"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_severe
  
  
  #AZ Mild
  
  filter(Andrews_Data, label == "Delta_Mild_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower", "#6 - waned upper",
                          "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Mild_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - AZ - Mild Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease",
               "#6 - Waned increased",
               "#7 - Waned decreased"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.3, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_mild
  
  #PF Mild
  
  filter(Andrews_Data, label == "Delta_Mild_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower", "#6 - waned upper",
                          "#7 - waned lower")
  d_simulated <- filter(Sim_Data, label == "Delta_Mild_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Delta - PF - Mild Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e",
                                  "#d8b365",
                                  "#c7eae5"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease",
               "#6 - Waned increased",
               "#7 - Waned decreased"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.3, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_mild
  
  
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







plot_f6b <- function(f4_data_alpha){
  
  #First, provide the data being fit to:
  
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
  
  
  #sensitivity #4
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
                         outcome = rep("mild", 4),
                         vaccine = rep("AZ",4),
                         variant = rep("Alpha", 4),
                         median = c(mean(waned_VE[(14-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(202-13)])),
                         lower = rep(NA, 4),
                         upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #5
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
                         outcome = rep("mild", 4),
                         vaccine = rep("AZ",4),
                         variant = rep("Alpha", 4),
                         median = c(mean(waned_VE[(14-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(202-13)])),
                         lower = rep(NA, 4),
                         upper = rep(NA, 4))
  
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #2
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
                         outcome = rep("mild", 4),
                         vaccine = rep("AZ",4),
                         variant = rep("Alpha", 4),
                         median = c(mean(waned_VE[(14-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(202-13)])),
                         lower = rep(NA, 4),
                         upper = rep(NA, 4))
  
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  #sensitivity #3
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
                         outcome = rep("mild", 4),
                         vaccine = rep("AZ",4),
                         variant = rep("Alpha", 4),
                         median = c(mean(waned_VE[(14-13):(69-13)]),
                                    mean(waned_VE[(70-13):(104-13)]),
                                    mean(waned_VE[(105-13):(139-13)]),
                                    mean(waned_VE[(140-13):(202-13)])),
                         lower = rep(NA, 4),
                         upper = rep(NA, 4))
  
  Sim_Data <- rbind(Sim_Data, Sim_Data_hold)
  
  
  
  
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
  
  #sensitivity #4
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
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
  
  #sensitivity #5
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
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
  
  #sensitivity #2
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
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
  
  #sensitivity #3
  waned_VE <- mild_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Mild_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
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
  
  
  
  
  
  
  ###############
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
  
  
  
  
  #sensitivity #4
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
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
  
  #sensitivity #5
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
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
  
  #sensitivity #2
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
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
  
  #sensitivity #3
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
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
  
  
  
  
  #sensitivity #4
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
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
  
  #sensitivity #5
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
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
  
  #sensitivity #2
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
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
  
  #sensitivity #3
  waned_VE <- severe_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Severe_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
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
  
  
  
  
  
  
  
  #sensitivity #4
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
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
  
  #sensitivity #5
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
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
  
  #sensitivity #2
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
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
  
  #sensitivity #3
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)+4])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_PF", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
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
  
  
  
  
  
  
  
  #sensitivity #4
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_4_all_upper[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#4 - all upper", 4),
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
  
  #sensitivity #5
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_5_all_lower[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#5 - all lower", 4),
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
  
  #sensitivity #2
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_2_mild_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#2 - reduced mild", 4),
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
  
  #sensitivity #3
  waned_VE <- death_plotting(c(24, as.numeric(f4_data_alpha$sens_3_all_absolute_reduction[c(2,10,18,3,11,19)])))
  
  Sim_Data_hold <- data.frame(label = rep("Alpha_Death_AZ", 4),
                         timeframe = c("2-9", #Days 14-69
                                       "10-14", #Days 70-104
                                       "15-19", #Days 105-139
                                       "20+"#Days 140-202
                         ),
                         sensitivity = rep("#3 - 0.1 reduction all", 4),
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
  
  
  
  Sim_Data$sensitivity <- factor(Sim_Data$sensitivity,
                                 levels = c("baseline",
                                            "#3 - 0.1 reduction all", "#4 - all upper",        
                                            "#5 - all lower"   ))
  
  
  ############################################
  ## PLOTTING
  ############################################
  
  
  
  #AZ Death
  
  filter(Andrews_Data, label == "Alpha_Death_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Death_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - AZ - Death") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    # scale_color_lancet(labels = c("#1 - Baseline",
    #                               "#3 - Absolute decrease",
    #                               "#4 - Relative increase",
    #                               "#5 - Relative decrease"
    #                               )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_death
  
  #PF Death
  
  filter(Andrews_Data, label == "Alpha_Death_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Death_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - PF - Death") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_death
  
  #AZ Severe
  
  filter(Andrews_Data, label == "Alpha_Severe_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Severe_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  data_hold <- rbind(d_Andrews, d_simulated)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - AZ - Severe Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_severe
  
  #PF Severe
  
  filter(Andrews_Data, label == "Alpha_Severe_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- sensitivities_keep <- c("baseline", "#4 - all upper",
                                                "#3 - 0.1 reduction all",
                                                "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Severe_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - PF - Severe Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.6, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_severe
  
  
  #AZ Mild
  
  filter(Andrews_Data, label == "Alpha_Mild_AZ") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Mild_AZ")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - AZ - Mild Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.3, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> AZ_mild
  
  #PF Mild
  
  filter(Andrews_Data, label == "Alpha_Mild_PF") %>%
    mutate(timeframe = fct_reorder(timeframe, desc(median))) -> d_Andrews
  
  sensitivities_keep <- c("baseline", "#4 - all upper", 
                          "#3 - 0.1 reduction all",
                          "#5 - all lower")
  d_simulated <- filter(Sim_Data, label == "Alpha_Mild_PF")
  d_simulated <- filter(d_simulated, sensitivity %in% sensitivities_keep)
  
  ggplot(d_simulated) +
    geom_point(data = d_Andrews, aes(x = timeframe, y = median,
                                     shape = sensitivity, group = sensitivity, shape = sensitivity),
               stroke = 1.5,
               color = "#841744",
               size = 2,
               key_glyph = "pointrange") +
    geom_errorbar(data = d_Andrews,
                  aes(x= timeframe, ymin = lower, ymax = upper,
                      group = sensitivity),
                  color = "#841744",
                  width = 0.3) +
    geom_line(aes(x = timeframe, y = median, color = sensitivity,
                  group = sensitivity),
              alpha = 0.9, size = 1) +
    theme_classic() +
    xlab("Weeks since second dose") + 
    ggtitle("Alpha - PF - Mild Disease") +
    ylab("VE") +
    labs(colour="Sensitivity") +
    scale_color_manual(values = c("darkgray", 
                                  "#5ab4ac",
                                  "#8c510a",
                                  "#01665e"
    ),
    labels = c("#1 - Baseline",
               "#3 - Absolute decrease",
               "#4 - Relative increase",
               "#5 - Relative decrease"
    )) +
    theme(axis.text=element_text(size=rel(1.2)),
          axis.title=element_text(size=rel(1.3)),
          axis.text.x = element_text(angle = 30, hjust=0.9),
          legend.text = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.3)),
          legend.spacing.y = unit(0.5, 'lines'),
          strip.text = element_text(size=rel(1.3))) +
    scale_shape_manual(name = "",
                       labels = c("Data"),
                       values= c(16)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3),
                       limits = c(0.3, 1)) +
    guides(color = guide_legend(byrow = TRUE,
                                order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 0.7))) -> PF_mild
  
  
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


