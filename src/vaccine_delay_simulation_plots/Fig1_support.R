prep_f1p1_data <- function(Real_vaccination, Counterfactual_vaccination) {
  
  #Prep the data frame
  cum_doses_by_age <- data.frame( total_doses= NA, date = NA, age_group = NA, dose = NA, category = NA)
  
  age_group_names <- c('0-4', '5-9', '10-14', '15-19',
                       '20-24', '25-29', '30-34', '35-39',
                       '40-44', '45-49', '50-54', '55-59',
                       '60-64', '65-69', '70-74', '75-79',
                       '80+', 'CHW', 'CHR')
  regions <- names(combined$base)
  
  #Fill the data frame
  for(i in 1:length(regions)){ #for each region
    for(j in 1:length(Real_vaccination[[i]]$doses[,1,1])){ #for each age group
      hold_real_doses1 <- data.frame(total_doses = (cumsum(Real_vaccination[[i]]$doses[j,1,])),
                                     date = seq(from = Real_vaccination[[i]]$date,
                                                length.out = length(Real_vaccination[[i]]$doses[1,1,])),
                                     age_group = rep(age_group_names[j],length(Real_vaccination[[i]]$doses[1,1,]) ),
                                     dose = rep(1, length(Real_vaccination[[i]]$doses[1,1,])),
                                     category = rep('Real', length(Real_vaccination[[i]]$doses[1,1,])))
      
      hold_real_doses2 <- data.frame(total_doses = (cumsum(Real_vaccination[[i]]$doses[j,2,])),
                                     date = seq(from = Real_vaccination[[i]]$date,
                                                length.out = length(Real_vaccination[[i]]$doses[1,2,])),
                                     age_group = rep(age_group_names[j],length(Real_vaccination[[i]]$doses[1,2,]) ),
                                     dose = rep(2, length(Real_vaccination[[i]]$doses[1,2,])),
                                     category = rep('Real', length(Real_vaccination[[i]]$doses[1,2,])))
      
      hold_real_doses_both <- hold_real_doses2
      hold_real_doses_both$total_doses <- hold_real_doses1$total_doses + hold_real_doses2$total_doses
      hold_real_doses_both$dose <- rep('All', length(hold_real_doses2$dose))
      hold_real_doses_both$category <- rep('All', length(hold_real_doses2$category))
      
      hold_doses <- rbind(hold_real_doses1, hold_real_doses2)
      hold_doses <- rbind(hold_doses, hold_real_doses_both)
      
      hold_counterfactual_doses1 <- data.frame(total_doses = (cumsum(Counterfactual_vaccination[[i]]$doses[j,1,])),
                                               date = seq(from = Counterfactual_vaccination[[i]]$date,
                                                          length.out = length(Counterfactual_vaccination[[i]]$doses[1,1,])),
                                               age_group = rep(age_group_names[j],length(Counterfactual_vaccination[[i]]$doses[1,1,]) ),
                                               dose = rep(1, length(Counterfactual_vaccination[[i]]$doses[1,1,])),
                                               category = rep('Counterfactual', length(Counterfactual_vaccination[[i]]$doses[1,1,])))
      
      hold_counterfactual_doses2 <- data.frame(total_doses = (cumsum(Counterfactual_vaccination[[i]]$doses[j,2,])),
                                               date = seq(from = Counterfactual_vaccination[[i]]$date,
                                                          length.out = length(Counterfactual_vaccination[[i]]$doses[1,2,])),
                                               age_group = rep(age_group_names[j],length(Counterfactual_vaccination[[i]]$doses[1,2,]) ),
                                               dose = rep(2, length(Counterfactual_vaccination[[i]]$doses[1,2,])),
                                               category = rep('Counterfactual', length(Counterfactual_vaccination[[i]]$doses[1,2,])))
      
      hold_counterfactuals <- rbind(hold_counterfactual_doses1, hold_counterfactual_doses2)
      hold_doses <- rbind(hold_doses, hold_counterfactuals)
      
      
      cum_doses_by_age <- rbind(cum_doses_by_age, hold_doses)
      
    }
  }
  
  #Remove the placeholder row
  cum_doses_by_age <- cum_doses_by_age[-1,]
  
  cum_doses_by_age$date <- sircovid::sircovid_date_as_date(cum_doses_by_age$date)
  cum_doses_by_age$age_group <- factor(cum_doses_by_age$age_group, levels = age_group_names)
  cum_doses_by_age$dose <- as.factor(cum_doses_by_age$dose)
  cum_doses_by_age$category <- as.factor(cum_doses_by_age$category)
  
  #This is now the full data set
  
  ####################
  #Now minimise the data, first sum over all regions
  
  cum_doses_by_age <- aggregate(cum_doses_by_age$total_doses, by=list(date=cum_doses_by_age$date,
                                                                      age_group = cum_doses_by_age$age_group,
                                                                      dose = cum_doses_by_age$dose,
                                                                      category = cum_doses_by_age$category), FUN=sum)
  colnames(cum_doses_by_age) <- c('date', 'age_group', 'dose', 'category', 'total_doses')
  
  #######
  #Next, sum over all ages: 
  cum_doses <- aggregate(cum_doses_by_age$total_doses, by=list(date=cum_doses_by_age$date,
                                                               dose = cum_doses_by_age$dose,
                                                               category = cum_doses_by_age$category), FUN=sum)
  colnames(cum_doses) <- c('date', 'dose', 'category', 'total_doses')
  
  #Now cut out the first doses:
  cum_doses <- filter(cum_doses, dose != '1')
  
  return(cum_doses)
  
}

 prep_f1p2_data <- function(simulation_output, fit_trajectories, combined){
   simulation_state <- simulation_output$state
   simulation_state$quantile <- as.character(simulation_state$quantile)
   fit_state <- fit_trajectories$state
   fit_state <- filter(fit_state, region == 'england')
   keep_quantiles <- c('2.5%', '50%', '97.5%')
   simulation_state <- filter(simulation_state, quantile %in% keep_quantiles)
   
   simulation_state <- filter( simulation_state, region == 'england')
   simulation_state <- subset(simulation_state, select=-c(vaccine_daily_doses,
                                                          booster_daily_doses,
                                                          beta_step,
                                                          scenario, group,
                                                          analysis, vaccine_status))
   fit_state <- subset(fit_state, select=-c(group, vaccine_status))
   simulation_state$type <- 'counterfactual'
   
   All_trajectories <- rbind(fit_state, simulation_state)
   All_trajectories <- filter(All_trajectories, date < as.Date('2021-09-14'))
   All_trajectories <- filter(All_trajectories, date > as.Date('2020-12-14'))
   #also: 
   keep_states <- c('diagnoses_admitted', 'diagnoses_admitted_inc', 'deaths', 'infections')
   All_trajectories <- filter(All_trajectories, state %in% keep_states) #Set state of interest here
   All_trajectories$quantile <- as.factor(All_trajectories$quantile)
   All_trajectories <- pivot_wider(All_trajectories, names_from = quantile, values_from = value)

   
   #We want the counterfactual admissions_cumulative to start from the same point as the fits
   for(i in c("diagnoses_admitted", "deaths", "infections")){
   #50%
   All_trajectories$`50%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] <- All_trajectories$`50%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] - All_trajectories$`50%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]  
   All_trajectories$`50%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] <- All_trajectories$`50%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] - All_trajectories$`50%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]
   #2.5%
   All_trajectories$`2.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] <- All_trajectories$`2.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] - All_trajectories$`2.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]    
   All_trajectories$`2.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] <- All_trajectories$`2.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] - All_trajectories$`2.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]
   #97.5%
   All_trajectories$`97.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] <- All_trajectories$`97.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i )] - All_trajectories$`97.5%`[which(All_trajectories$type == 'counterfactual' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]
   All_trajectories$`97.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] <- All_trajectories$`97.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i )] - All_trajectories$`97.5%`[which(All_trajectories$type == 'fit' & All_trajectories$state == i & All_trajectories$date == as.Date("2020-12-15") )]
   }
   #Prep the real world data to plot over.
   Hosp_data <- data.frame(date = combined$data$east_of_england$date_end,
                           all_admissions = rep(0,547))
   for( i in 1:7){
     Hosp_data$all_admissions <- Hosp_data$all_admissions + combined$data[[i]]$all_admission
     #set NAs to 0:
     Hosp_data$all_admissions[is.na(Hosp_data$all_admissions)] = 0
   }
   Hosp_data$date <- as.Date(sircovid::sircovid_date_as_date(Hosp_data$date))
   #Remove outside of time window we're interested in
   Hosp_data <- filter(Hosp_data, date > as.Date('2020-12-15'))
   
   
   #We can output the cumulative data here to save computation.
   #This is some over-writing to make plotting f1p3 easier -
   # Ribbon_lower is NOT the lower of the CI, it's for marking a stacked plot where factual and counterfactual meet.
   Cumulative_data <- All_trajectories
   Cumulative_data <- filter(Cumulative_data, date > as.Date('2020-12-15'))
   names(Cumulative_data)[names(Cumulative_data) == "2.5%"] <- 'Ribbon_lower'
   Cumulative_data$Ribbon_lower[which(Cumulative_data$state == 'diagnoses_admitted' & Cumulative_data$type == 'counterfactual')] <- Cumulative_data$`50%`[which(Cumulative_data$state == 'diagnoses_admitted' & Cumulative_data$type == 'fit')]
   Cumulative_data$Ribbon_lower[which(Cumulative_data$state == 'diagnoses_admitted' & Cumulative_data$type == 'fit')] <- 0
   
   ret <- list("trajectories" = All_trajectories, "Real_data" = Hosp_data, "Cumulative_data" = Cumulative_data)
   
   return(ret)
   
 }
 
 
 
 
 
 ##############
 #Plotting functions
 
 plot_f1p1 <- function(f1p1_data){
   
   #Key dates here:
   grey_lines <- c(
     "2021-01-05", ## 16. Lockdown 3 starts
     "2021-03-08", ## 17. Step 1 of roadmap: schools reopen
     "2021-04-19", ## 19. Step 2 of roadmap: outdoors hospitality (04-12) 
     ##     and schools return (04-19)
     "2021-05-17", ## 20. Step 3 of roadmap: indoors hospitality
     "2021-07-19") ## 24. Step 4
   grey_labels <- as.Date(grey_lines)
   #We give some label positions a nudge to aid plotting
   grey_labels[5] <- grey_labels[5] 
   grey_labels[4] <- grey_labels[4] + 10
   
   ggplot(data = f1p1_data, aes(x = date, y= total_doses/1000000, color = category)) +
     geom_line(size = 1.25) + theme_classic() + ylab('Cumulative vaccines \nadministered (millions)') +
     geom_ribbon(data = filter(f1p1_data, category == 'All'), aes(ymin = 0, ymax = total_doses/1000000 , fill = category)) +
     geom_ribbon(data = filter(f1p1_data, category == 'Counterfactual'), aes(ymin = 0, ymax = total_doses/1000000 , fill = category)) +
     geom_ribbon(data = filter(f1p1_data, category == 'Real'), aes(ymin = 0, ymax = total_doses/1000000 , fill = category)) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
     scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
     scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
     guides(fill = guide_legend(byrow = TRUE)) +
     labs(fill = "Vaccine doses") + theme(legend.spacing.y = unit(0.5, 'cm')) +
     coord_cartesian(xlim = c(as.Date('2020-12-01'),as.Date('2021-09-15'))) +
     scale_fill_discrete(labels = c('All',
                                    'Dose 2 \n(Counterfactual)',
                                    'Dose 2 \n(Model fit)')) +
     geom_vline(xintercept = as.Date('2021-03-31'), alpha = 0.9, color = 'black') +
     geom_vline(xintercept = as.Date(grey_lines), alpha = 0.7, color = 'gray23', lty = 'dashed') +
     annotate('text', x=as.Date('2021-03-25'), label="Delta emergence", y=50, colour="black", size = 5, angle = 90) +
     scale_color_manual(values = c("gray80", "#9BC362", "#B25A82"), guide = "none") +
     scale_fill_manual(breaks = c("All", "Real", "Counterfactual"),
                       values = c("gray80", "#B25A82", "#9BC362"),
                       labels = c('All',
                                  'Dose 2 \n(12-week strategy)',
                                  'Dose 2 \n(3-week strategy)'
                                  ))+ xlab(NULL) +
      geom_label(data = data.frame(dates = grey_labels, ypos = rep(75,5)),
                 aes(x = dates, y = ypos, label = c('Lockdown 3 starts', 'Roadmap \nStep 1',
                                                    'Roadmap \nStep 2', 'Roadmap \nStep 3',
                                                    'Roadmap \nStep 4'), color = NULL), size = 3 ) +
      theme(axis.text=element_text(size=rel(1.2)),
           axis.title=element_text(size=rel(1.3)),
           legend.text = element_text(size=rel(1.2)),
           legend.title = element_text(size=rel(1.3))) -> p1
   
   return(p1)
   
 }
 
 
 plot_f1p2 <- function(f1p2_data){
   
   #Key dates here:
   grey_lines <- c(
     "2021-01-05", ## 16. Lockdown 3 starts
     "2021-03-08", ## 17. Step 1 of roadmap: schools reopen
     "2021-04-19", ## 19. Step 2 of roadmap: outdoors hospitality (04-12) 
     ##     and schools return (04-19)
     "2021-05-17", ## 20. Step 3 of roadmap: indoors hospitality
     "2021-07-19") ## 24. Step 4
   
   ggplot(data = filter(f1p2_data$trajectories, state == 'diagnoses_admitted_inc'), aes(x= date, y = `50%` )) +
     geom_line(size = 1, aes(color = type)) + 
     geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = type), alpha = 0.3, color = NA, show.legend = FALSE) + 
     theme_classic() + ylab('Daily hospital admissions') +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
     
     scale_fill_manual(values=c("#9BC362", "#B25A82")) +
      scale_color_manual(breaks = c("fit", "counterfactual"),
                        values = c("#B25A82", "#9BC362"),
                        labels = c("12-week", '3-week')) +
      coord_cartesian(xlim = c(as.Date('2020-12-01'),as.Date('2021-09-15'))) +
     labs(color = "Strategy") +
     geom_vline(xintercept = as.Date('2021-03-31'), alpha = 0.9, color = 'black') +
     geom_vline(xintercept = as.Date(grey_lines), alpha = 0.7, color = 'gray23', lty = 'dashed') +
     xlab(NULL) +
     theme(axis.text=element_text(size=rel(1.2)),
           axis.title=element_text(size=rel(1.3)),
           legend.text = element_text(size=rel(1.2)),
           legend.title = element_text(size=rel(1.3))) +

     geom_point(data = f1p2_data$Real_data, aes(x = date, y = all_admissions), alpha = 0.7, shape = 18,
                color = 'black')  -> p2
   
   return(p2)
 }
 
 plot_f1p3 <- function(f1p3_data){
   
    #Key dates here:
    grey_lines <- c(
       "2021-01-05", ## 16. Lockdown 3 starts
       "2021-03-08", ## 17. Step 1 of roadmap: schools reopen
       "2021-04-19", ## 19. Step 2 of roadmap: outdoors hospitality (04-12) 
       ##     and schools return (04-19)
       "2021-05-17", ## 20. Step 3 of roadmap: indoors hospitality
       "2021-07-19") ## 24. Step 4
    grey_labels <- as.Date(grey_lines)
    #We give some dates a nudge to aid plotting
    grey_labels[5] <- grey_labels[5] 
    grey_labels[4] <- grey_labels[4] + 10
    
    ggplot() +
       #cumulative
       geom_ribbon(data = filter(f1p3_data, state == 'diagnoses_admitted' & type == 'fit'),
                   aes(x = date, ymin = Ribbon_lower/1000, ymax = `50%`/1000 , fill = type )) +
       geom_ribbon(data = filter(f1p3_data, state == 'diagnoses_admitted' & type == 'counterfactual'),
                   aes(x = date, ymin = Ribbon_lower/1000, ymax = `50%`/1000 , fill = type )) +
                      theme_classic() + ylab('Cumulative hospital \nadmissions (thousands)') +
       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
       scale_x_date(date_breaks = "1 month", date_labels = "%b %y", name = 'Date') +
       scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
       scale_fill_manual(breaks = c("fit", "counterfactual"),
                         values = c("#B25A82", "#9BC362"),
                         labels = c("12-week", '3-week')) +
       coord_cartesian(xlim = c(as.Date('2020-12-01'),as.Date('2021-09-15'))) +
      
       labs(fill = "Strategy") +
       geom_vline(xintercept = as.Date('2021-03-31'), alpha = 0.9, color = 'black') +
       geom_vline(xintercept = as.Date(grey_lines), alpha = 0.7, color = 'gray23', lty = 'dashed') +
       theme(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.text = element_text(size=rel(1.2)),
             legend.title = element_text(size=rel(1.3)))  -> p3
    
    return(p3)
    
    
 }
 
 