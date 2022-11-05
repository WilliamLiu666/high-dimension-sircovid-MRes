
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------

#------------------------------------------------------------------------------
#DATA SETUP TASK
fits_data <- orderly::orderly_run(
  "vaccine_delay_fits_data"
)

orderly::orderly_commit(fits_data)
#------------------------------------------------------------

## Parameters fits task sets up before fitting for each region
pars_bb <- orderly::orderly_run(
  "vaccine_delay_parameters_fits",
  parameters = list(deterministic = FALSE),
  use_draft = "newer")
orderly::orderly_commit(pars_bb)


### ---------------------------------------------------------------------------
##--------------------
## SHORT RUNS
# Set as TRUE for short runs to investigate and de-bug
# Set as FALSE for far longer runs and better fits
short_run <- FALSE

##DETERMINISTIC
# Model can be run deterministically for quicker de-bugging process or fit tuning
# Note that simulation tasks downstream require stochastic (deterministic <- FALSE)
deterministic <- FALSE
##--------------------
regions <- sircovid::regions("england")

#-------------
#RUN
bb_fits_ids <- lapply(X = regions,
                 FUN = function(x) {
                   orderly::orderly_run('vaccine_delay_fits',
                                        parameters = list(region = x,
                                                          short_run = short_run,
                                                          deterministic = deterministic),
                                        use_draft = "newer")})

lapply(X= bb_fits_ids,
       FUN = function(x){
         orderly::orderly_commit(x)
       })


# combine regions
combined_bb <- orderly::orderly_run('vaccine_delay_fits_combined',
                                    parameters = list(short_run = short_run,
                                                      deterministic = deterministic
                                                      ),
                                    use_draft = "newer")

orderly::orderly_commit(combined_bb)




##############################################################
#SIMULATION PIPELINE
################################################################
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------
########################
#Sensitivity Parameters

#These parameters can be used to investigate specific sensitivity cases
#as displayed in the supplementary material.

#Changing the VE assumptions.
#FALSE is the default, with no differences to VE in the counterfactual
#"first_dose" alters only first dose VE as specified in sensitivity 9 in Table S9
#"targetted" alters specifically second dose and reduced protection VE (sensitivities 2/3 in Table S9)
# Setting a specific value will add this value to all counterfactual VEs
vaccine_efficacy_alteration <- FALSE
#vaccine_efficacy_alteration <- "first_dose"
#vaccine_efficacy_alteration <- "targetted"
#vaccine_efficacy_alteration <- -0.1

#Change the time to wane
#Default is FALSE, waning is kept at an average of 24 weeks
#Set to a specific value (in days) for custom mean waning time
mean_waning_time <- FALSE

#Change time between doses (in days)
#Default for all counterfactuals is 3 weeks:
days_between_doses <- 3*7

########################

## Set fit data range
data_date <- "2020-12-08"
end_date <- "2022-05-01" #Fits use a longer data period

#Number of particles to use in counterfactual simulation. Greater value increases runtime.
n_par <- 200
#n_par <- 25


param_simulation <- "vaccine_delay_simulation_parameters"
simulation <- "vaccine_delay_simulation"
# ---------------------------------

## 1. Simulation parameters ----
# Run
sim_params <- orderly::orderly_run(param_simulation,
                     parameters = list(date_of_fits = data_date,
                                       simulation_end = end_date),
                     use_draft = "newer")

orderly::orderly_commit(sim_params)
# ----

## 2. Simulation task ----

# Run
sim_task <- orderly::orderly_run(simulation,
                     parameters = list(date_of_fits = data_date,
                                       n_par = n_par,
                                       vaccine_efficacy_alteration = vaccine_efficacy_alteration,
                                       days_between_doses = days_between_doses, 
                                       mean_waning_time = mean_waning_time,
                                       short_run = short_run,
                                       debug = TRUE),
                     use_draft = "newer")

orderly::orderly_commit(sim_task)

###########
#PLOTS

sim_plots <- orderly::orderly_run('vaccine_delay_simulation_plots',
                     parameters = list(date_of_fits = data_date,
                                       days_between_doses = days_between_doses,
                                       vaccine_efficacy_alteration = vaccine_efficacy_alteration,
                                       n_par = n_par,
                                       short_run = short_run),
                     use_draft = "newer")

orderly::orderly_commit(sim_plots)


