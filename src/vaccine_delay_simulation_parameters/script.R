source("global_util.R")
source("global_vaccine.R")

version_check("sircovid", "0.13.15")
version_check("spimalot", "0.7.11")

date <- as.Date(date_of_fits)

## Assumptions that never change
# Seasonality (sine wave)
rt_seasonality_date_peak <- sircovid::sircovid_date("2020-02-15")
rt_seasonality_amplitude <- 0
# Log-normal standard deviation for generating Rt; 
# this will only be used if not specified in the input csv file
# TO DO: change this so it's clearer. Maybe remove default value from here?
rt_sd <- 0.0001
# Multiplicative amount Rt decreases when schools close (0.15 is a 15%
# reduction, so a 0.85 multiplier to contacts with school open).
# Applied multiplicatively with seasonality
rt_schools_modifier <- 0
# Schools open/close schedule
rt_schools_schedule <- read_csv("schools_schedule.csv")
rt_scenarios <- read_rt_scenarios(file.path("mtp", "rt_scenarios.csv"),
                                  date, rt_sd)

future_doses <- read.csv(file.path("mtp", "vaccine_schedule.csv"),
                         stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")

## This parameter (vaccine_eligibility_min_age)
# set below is used in both ‘mtp’ and ‘add omicron’ simulations.
vaccine_eligibility_min_age <- 5

simulation_regions <- sircovid::regions("england")

## Ensure that the MTP simulations always start on a Monday
output_day_start <- "Mon"
week <- as.Date(date) + 0:6
date_simulation_start <-
  week[lubridate::wday(week, label = TRUE) == output_day_start]

if (is.numeric(simulation_end)) {
  date_simulation_end <- as.Date(date_simulation_start) + simulation_end
} else {
  date_simulation_end <- as.Date(simulation_end)
}

future_doses <- FIXME_VACCINE_SCHEDULE(future_doses, simulation_regions)

# A vector of names of simulation outputs to keep from the
# simulation; this list should be kept fairly small.
output_keep <- c("deaths", "deaths_inc", "admitted", "diagnoses",
                  "infections", "infections_inc", "hosp", "react_pos")

grid <- tidyr::expand_grid(
  vaccine_daily_doses = "central",
  booster_daily_doses = "central",
  beta_step = unique(rt_scenarios$scenario))

## Not totally sure that either of these are useful, but here we
## are.
grid$scenario <- grid$beta_step
grid$analysis <- sprintf("analysis %d", seq_len(nrow(grid)))

# This basically combines the parameters we change for the simulation
parameters <- list(
  vaccine_daily_doses = future_doses$vaccine_daily_doses,
  booster_daily_doses = future_doses$booster_daily_doses,
  vaccine_eligibility_min_age = vaccine_eligibility_min_age,

  ## Simulation assumptions and parameters
  rt_sd = rt_sd,
  rt_schools_schedule = rt_schools_schedule,
  rt_schools_modifier = rt_schools_modifier,
  rt_scenarios = rt_scenarios,
  rt_seasonality_date_peak = rt_seasonality_date_peak,
  rt_seasonality_amplitude = rt_seasonality_amplitude)

## Here, write out the list of parameters that will be provided.
## This wil match the names(parameters) but it is best to keep this
## separate as a check that you have provided the full list and that
## everyone can see what has been added.
expected <- c("vaccine_eligibility_min_age",
              "vaccine_daily_doses",
              "booster_daily_doses",
              "beta_step",
              "rt_sd",
              "rt_schools_schedule",
              "rt_schools_modifier",
              "rt_scenarios",
              "rt_seasonality_date_peak",
              "rt_seasonality_amplitude")

output <- spimalot::spim_simulate_control_output(
  keep = output_keep, time_series = TRUE,
  rt = TRUE, rt_weighted = TRUE, rt_type = c("Rt_general", "eff_Rt_general"),
  vaccination = TRUE,
  #NEW STATE BY AGE (if FALSE won't run)
  state_by_age = TRUE) 

spim_state_names <- read_csv("spim_state_names.csv")
spim_state_names <- setNames(spim_state_names$spim_name,
                              spim_state_names$sircovid_name)

control <- spimalot::spim_simulate_control(
  flavour = "mtp", regions = simulation_regions,
  date_start = date_simulation_start, date_end = date_simulation_end,
  expected = expected, parameters = parameters, grid = grid, output = output)

saveRDS(control, "spim_parameters_simulate.rds")
