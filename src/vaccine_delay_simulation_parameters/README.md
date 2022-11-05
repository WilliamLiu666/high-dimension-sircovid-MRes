### ---- README - Documentation to add or modify simulation scenarios ---- ###

** THIS IS OUT OF DATE AND NEEDS REWRITING **


Producing simulation scenarios for MTPs of CO asks requires updating:
1. Corresponding section in `script.R`
2. Vaccine doses available in `mtp_vaccine_schedule.csv` and/or
   `co_vaccine_schedule.csv`
3. If additioanl school closures are to be modelled (e.g. lockdown), the
   schedule will need modifying in `schools_schedule.csv`.
   
This documentation will focus on updating specific sections in `script.R`.
For the case of the vaccine schedule files consult SPI-M email (usually
circulated by Marc Baguelin at COP on a Friday). Keep in mind we can only
simulate doses that will be made available on or after the date the simulation
will be started. Often asks include doses in the past; retrospective allocation
of doeses is not supported in our model framework. For the schools schedule,
please refer directly to the specified file. It is currently up to date with the
regular school calendar for all nations until July 2022.


### I. Set-up simulations parameters for MTPs ----

The task `rtm_parameters_simulate` requires specifying orderly parameters:
- date_of_fits: character string for the date used for running
  rtm_spim_restart_fits_combined
- assumptions: either of "central", "pessimistic" or "optimistic"
- model_type: either of "BB" or "NB"
- simulation_end: can either be numeric or character string. IF numeric, this
  will be the number of days to simulate for, after date_of_fits. If character
  string, this should be a specific date in which the simulation will end

When `simulatew what == "mtp"`, first section of script.R will need updating

## 1. Basic variables

# a. Basic scenario variables
This is the only section that, for most cases, will need updating.
- scenario_names: SPI-M is particular about the naming of scenarios. Check
  weekly email and ask for clarification if needed (Matt and Alaistair are usual
  points of call in SPI-M secretariat). Default to "MTP" if the ask is for
  `normal MTPs only`; this is a basic scenario where Rt (beta) is assumed to be
  constant during the simumaltion. Else, specify a concatenated vector
  corresponding to scenarios (e.g. c("MTP 0.9", "MTP 1.1", "MTP 1.3"))
- n_gradual_steps: integer of the number of steps for the scenario. If only 
  simulating "MTP normal", this should be one; if simulating scenarios with one
  change, this is two, etc.
- date_changein_rt: a character string or vector of strings for the dates at
  which Rt will change 
- scaling: logical, specify whether Rt change will be a scaling (up/down) of
  current Rt. If FALSE, scenarios will correspond to actual target Rt values

# b. Regions to simulate
For MTPs, this should be all. Currently set to "england", as simulating for the
whole of the UK remains broken. Issue raised in sircovid (#365).

# c. Simulation start and end dates
99.9% of MTPs start of the Monday of the week of the ask, hence
date_simulation_start will seek the next Monday on or after the date of last
data point.


## 2. Future doses schedules

This section reads-in "mtp_vaccine_schedule.csv" and uses regional_pop to break
down doses available by region over time. Only part that needs editing is the
actual CSV, using schedules suplied by SPI-M. See note above.


## 3. Define simulation scenarios

# a. Rt scenarios
This section will also need frequent updating, depending on the ask. Bellow are
a number of examples for the most frequent styles of asks:

- MTP normal - again, this is a scenario with no change in Rt, rather a constat
  Rt that accounts for a future doses schedule, schools schedule and seasonality
  
  ## where scenario_names <- "MTP" and scaling <- TRUE
  rt_scenarios <- list(
    rt_values = 1,
    date_rt_start = date_simulation_start,
    run_grid = tolower(scenario_names),
    scaling = scaling
  )
  
- Rt scenarios - a number of scenarios have been requested, setting Rt to 
  arbitrary values on a given date
  
  ## e.g. where a specific Rt values are requested
  ## scenario_names <- c("MTP 0.9", "MTP 1.1", "MTP 1.3", "MTP 1.4")
  ## scaling <- FALSE
  ## n_gradual_steps <- 2
  ## date_change_in_rt <- as.Date("2021-12-06")
  rt_scenarios <- list(
    rt_values = as.numeric(gsub("[^0-9.]", "",  scenario_names)),
    date_rt_start = as.Date("2021-12-06"),
    run_grid = tolower(gsub("[ .]", "_",  scenario_names)),
    scaling = scaling
  )
  
  
  ## e.g. where scaling changes in Rt at the start of the simulation requested
  ## scenario_names <- c("MTP 1", "MTP 1.1", "MTP 1.3")
  ## scaling <- TRUE
  ## n_gradual_steps <- 1
  ## date_change_in_rt <- date_simulation_start
  rt_scenarios <- list(
    rt_values = 1,
    date_rt_start = date_simulation_start,
    run_grid = tolower(scenario_names),
    scaling = scaling
  )
