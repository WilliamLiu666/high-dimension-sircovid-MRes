source("global_util.R")
source("global_vaccine.R")

version_check("sircovid", "0.13.15")
version_check("spimalot", "0.7.11")


control <- readRDS("spim_parameters_simulate.rds")

n_threads <- spimalot::spim_control_cores()
message(sprintf("Running on %d threads", n_threads))

date <- as.Date(date_of_fits)

## Load parameters and fit outputs
combined <- readRDS("combined.rds")
combined$steps_per_day <- 4 
combined$step <- sircovid::sircovid_date(combined$date) * combined$steps_per_day
combined$dt <- 1 / combined$steps_per_day

# Model all nations or just England, as defined in rtm_parameters_simulate
celtic_nations <- c("scotland", "northern_ireland", "wales")
if (all(celtic_nations %in% control$regions)) {
  nations <- sircovid::regions("nations")
} else {
  nations <- "england"
}

## Set seed_voc = FALSE for MTP but useful for future seeding
combined_processed <- spimalot::spim_simulate_prepare(
  combined, control$parameters, n_par,
  regions = control$regions,
  seed_voc = FALSE
)

#Output the existing vaccine schedule:
Real_vaccine_schedule <- list(east_of_england = combined_processed$base[[1]]$vaccine_schedule,
                              midlands = combined_processed$base[[2]]$vaccine_schedule,
                              london = combined_processed$base[[3]]$vaccine_schedule,
                              north_east_and_yorkshire = combined_processed$base[[4]]$vaccine_schedule,
                              north_west = combined_processed$base[[5]]$vaccine_schedule,
                              south_east = combined_processed$base[[6]]$vaccine_schedule,
                              south_west = combined_processed$base[[7]]$vaccine_schedule)

#We now overwrite the vaccine schedule in each region with a new one determined
#by the days_between_doses parameter
if(!isFALSE(days_between_doses)){
  message('Altering vaccine schedule to a counterfactual')
  #For every region...
  for(i in 1:length(control$regions)){
    #Total daily doses:
    daily_doses_value <- colSums(combined_processed$base[[i]]$vaccine_schedule$doses, dims = 1L)
    daily_doses_value <- colSums(daily_doses_value, dims = 1L)
    date_start <- combined_processed$base[[i]]$vaccine_schedule$date
    mean_days_between_doses <- days_between_doses
    
    vaccine_uptake <- combined_processed$base[[i]]$vaccine_uptake

    vaccine_uptake <- vaccine_uptake[, 1]
    vaccine_eligibility <-
      sircovid::vaccine_eligibility(combined_processed$base[[i]]$vaccine_eligibility_min_age)
    
    priority_population <- sircovid::vaccine_priority_population(region = control$regions[i], 
                                                                 uptake = vaccine_uptake * vaccine_eligibility)
    
    new_vaccine_schedule <- sircovid::vaccine_schedule_future(
      date_start,
      daily_doses_value,
      mean_days_between_doses,
      priority_population
    )
   #Then replace the existing vaccine schedule with a new one:
    combined_processed$base[[i]]$vaccine_schedule <- new_vaccine_schedule
    #and the effective schedule (factoring in time until) :
    new_vaccine_schedule_effect <-
      shift_doses(new_vaccine_schedule, combined_processed$base[[i]]$vaccine_days_to_effect)
    
    combined_processed$base[[i]]$vaccine_schedule_effect <- new_vaccine_schedule_effect
    
  }
}

Counterfactual_vaccine_schedule <- list(east_of_england = combined_processed$base[[1]]$vaccine_schedule,
                              midlands = combined_processed$base[[2]]$vaccine_schedule,
                              london = combined_processed$base[[3]]$vaccine_schedule,
                              north_east_and_yorkshire = combined_processed$base[[4]]$vaccine_schedule,
                              north_west = combined_processed$base[[5]]$vaccine_schedule,
                              south_east = combined_processed$base[[6]]$vaccine_schedule,
                              south_west = combined_processed$base[[7]]$vaccine_schedule)



rt_future <- control$parameters$rt_scenarios

write_csv(rt_future, "rt_future.csv")

### keep only parameters from the last epoch
combined_processed$pars[] <-
  lapply(combined_processed$pars, function(p) p[[length(p)]]$pars)
##The rel parameters now have dim 19x4x4

# Obtain population and % infected by region for post-processing
population <- spim_population(combined_processed, ignore_uk = TRUE)

prop_infected <- spimalot::spim_prop_infected(combined, population,
                                              regions = control$regions)


beta_step <- compute_beta_step(combined_processed, rt_future, control)
combined_processed$rt <- NULL
control <- spimalot::spim_simulate_set_beta_step(control, beta_step)

control$parameters$rt_seasonality_amplitude <- 0
control$parameters$rt_sd <- 0.001
control$parameters$rt_schools_modifier <- 0

## Logical simulation parameters, expanded over our grid of scenarios
grid_pars <- spimalot::spim_simulate_parameter_grid(control)

## A full set of nrow(grid) * n(regions) * n_pars parameters, each of
## which can create a sircovid model directly.
pars <- prepare_parameter_update(grid_pars, combined_processed, control)
rt_critical_dates <- rt_critical_dates(rt_future, combined)

#We wish to extract parameters relating to p_C and p_H
#p_C: 1 x 19(age)
p_C <- pars[[1]]$london[[1]]$p_C_step

#h_H:  7(region) x 1649(time) x 19(age) 
#h_H captures the piecewise linear change over time in hospitalisation risk conditional on symptomatic disease.

#Take the mean across all particles n_par
h_H <- array(0, dim = c(length(control$regions), #7
                        length(pars[[1]][[1]][[1]]$p_H_step[,1]), #1649
                        length(p_C))) #19

for(i in seq_len(length(control$regions))){
  for(j in seq_len(length(pars[[1]][[1]]))){  #n_par that gets actually used
    h_H[i,,] <- h_H[i,,] + pars[[1]][[i]][[j]]$p_H_step
  }
  h_H[i,,] <- h_H[i,,]/length(pars[[1]][[1]])  #Divided over total length of particles
}

#Then, rel_p_hosp_if_sympt is calculated for each strain from VE data, 
#and does not vary by region or particle

#rel_p_hosp_if_sympt_VE: (19x4(strain)x4(vacc))
rel_p_hosp_if_sympt_VE <- pars[[1]]$london[[1]]$rel_p_hosp_if_sympt

#Lastly, the strain_rel_p (1x4, strain), 
#which is different by region, and npar particle
#strain_rel_p_variant; (7(region) x 4 (strain))
strain_rel_p_variant <- array(0, dim = c(length(control$regions), #7
                                         length(rel_p_hosp_if_sympt_VE[1,,1]))) #4
for(i in seq_len(length(control$regions))){
  for(j in seq_len(length(pars[[1]][[1]]))){  #length of n_par that is used
    strain_rel_p_variant[i,] <- strain_rel_p_variant[i,] + pars[[1]][[i]][[j]]$strain_rel_p_hosp_if_sympt
  }
  strain_rel_p_variant[i,] <- strain_rel_p_variant[i,]/length(pars[[1]][[1]])  #divided over total number of particles
}

#Then put them all in one list.
p_C_and_H <- list(p_C = p_C,
                  h_H = h_H,
                  rel_p_hosp_if_sympt_VE = rel_p_hosp_if_sympt_VE,
                  strain_rel_p_variant = strain_rel_p_variant)

## Results of running the simulation over each of our scenarios.
seed <- 42
res <- spim_simulate_local(pars, combined_processed, rt_critical_dates,
                           control$date_end, control$output, n_threads,
                           seed)

dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)

#Output vaccine schedules:
saveRDS(Counterfactual_vaccine_schedule, "outputs/Counterfactual_vaccine_schedule.rds")
saveRDS(Real_vaccine_schedule, "outputs/Real_vaccine_schedule.rds")

#Output p_C and p_H
saveRDS(p_C_and_H, "outputs/p_C_and_H.rds")

message("Creating national trajectories")

res_england <- lapply(res, spimalot::spim_simulate_process_output,
                  combined_region = "england",
                  regions = sircovid::regions("england"),
                  incidence_states = "diagnoses_admitted",
                  reset_states = TRUE,
                  simulation_start_date = control$date_start)
res_england <- lapply(res_england, spimalot::spim_simulate_create_summary,
                  #at = seq(0.05, 0.95, 0.05))
                  at = c(0.025,0.05, 0.25, 0.5, 0.75, 0.95, 0.975))

summary_tidy <- spimalot::spim_simulate_tidy_states(res_england,
                                                    control$grid,
                                                    combined_processed)

# From here on, use summary_tidy (formatted output)
summary_tidy$population_infected <- list(population = population,
                                         prop_infected = prop_infected)
saveRDS(summary_tidy, "outputs/summary.rds")
####################
summary_restart_tidy <-
  spimalot::spim_simulate_process_output(
    combined$simulate, "england", control$regions,
    c("diagnoses_admitted", "deaths", "infections"),
    reset_states = FALSE, rm.rtUK = TRUE) %>%
  spimalot::spim_simulate_create_summary() %>%
  tidy_state_one_LOCAL(list(type = "fit"))


saveRDS(summary_restart_tidy, "outputs/summary_restart.rds")

summary_tidy_state <- summary_tidy$state %>%
  dplyr::filter(quantile %in% c("2.5%", "5%", "95%", "97.5%", "50%"))
plot_analysis <- summary_tidy_state$analysis[1]
plot_scenario <- summary_tidy_state$scenario[1]

summary_tidy_state_by_age <- summary_tidy$state_by_age 
saveRDS(summary_tidy_state_by_age, "outputs/summary_tidy_state_by_age.rds")
# For all uptake plots
doses <- spimalot::spim_calculate_doses(
  summary_tidy, population$england,
  plot_scenario
)

#We trim off everything after September 13th 2021 as this is the end point of fits.
summary_tidy_state <- filter(summary_tidy_state, date < as.Date('2021-09-14'))
summary_tidy_state_by_age <- filter(summary_tidy_state_by_age, date < as.Date('2021-09-14'))

## Plot Rt for strain 1, 2, and average. Plot to pdf to overflow to multiple
##  pages when needed
pdf("figs/check_Rt_all.pdf", width = 11, height = 7)
p_Rt_all <- spimalot::spim_plot_check_rt(summary_tidy_state,
                        unique(rt_future$date),
                        summary_restart_tidy$state,
                        #summary_tidy_state,
                        c("eff_Rt_general_strain_1",
                          "eff_Rt_general_strain_2",
                          "eff_Rt_general_both",
                          "Rt_general_strain_1",
                          "Rt_general_strain_2",
                          "Rt_general_both"))
plot(p_Rt_all)
dev.off()

## Plot average Rt. Plot to pdf to overflow to multiple
##  pages when needed
pdf("figs/check_Rt_avg.pdf", width = 11, height = 7)
p <- spimalot::spim_plot_check_rt(summary_tidy_state,
                                  unique(rt_future$date),
                                  summary_restart_tidy$state,
                                  #summary_tidy_state,
                                  c("eff_Rt_general_both",
                                    "Rt_general_both"))
plot(p)
dev.off()

## Plot deaths_inc, infections_inc, hosp, diagnoses_admitted_inc over time.
##  facet across states and analyses
pdf("figs/check_states_by_age.pdf", width = 11, height = 7)
p <- spimalot::spim_plot_check_state_by_age(summary_tidy_state_by_age,
                                      ana = 'analysis 1',
                                      scen = 'Fixed')
plot(p)
dev.off()

## Plot deaths_inc, infections_inc, hosp, diagnoses_admitted_inc over time.
##  facet across states and analyses
pdf("figs/check_states.pdf", width = 11, height = 7)
p_states <- spimalot::spim_plot_check_state(summary_tidy_state,
                                     summary_restart_tidy$state)
plot(p_states)
dev.off()

## Plot dose uptake by age group over time. Facet by analysis and dose.
pdf("figs/check_doses_by_age.pdf", width = 11, height = 7)
p <- spimalot::spim_plot_check_doses(doses)
plot(p+ xlim(as.Date('2020-12-14'), as.Date('2021-09-13')))
dev.off()

## Plot total daily doses over time. Facet by analysis.
pdf("figs/check_daily_doses.pdf", width = 11, height = 7)
p <- spimalot::spim_plot_check_total_doses(doses)
plot(p + xlim(as.Date('2020-12-14'), as.Date('2021-09-13')))
dev.off()

## Plot proportion of dose uptake over time and age. Facet by analysis and dose.
pdf("figs/check_proportion_doses.pdf", width = 11, height = 7)
p <- spimalot::spim_plot_check_uptake(doses)
plot(p+ xlim(as.Date('2020-12-14'), as.Date('2021-09-13')))
dev.off()

pdf("figs/check_daily_infections.pdf", width = 11, height = 7)
plot(1, 1, main = "simulate_flavour == 'mtp'")
dev.off()



