source("global_util.R")

version_check("sircovid", "0.14.7")
version_check("spimalot", "0.8.15")

date <- "2020-09-13"
model_type <- "BB"

#Some data streams are unreliable for the last few days. Define how many days here.
trim_deaths <- 4
trim_pillar2 <- 5

## MCMC control (only applies if short_run = FALSE)
burnin <- 500
n_mcmc <- 1500
chains <- 4
kernel_scaling <- 0.2

#Checks current region is valid
region <- spimalot::spim_check_region(region, FALSE)

#Checks and assembles all the outputs from vaccine_delay_parameter_fits task
#Output pars is a list containing:
#info - the loaded info.csv
#prior - the loaded prior.csv
#proposal - the loaded proposal.csv
#transform - the functions in transform.R
#raw - exact output of first 3 csvs, without any treatment
#base - the baseline.rds object
#mcmc - some sort of initialisation object built from the above to pass to the mcmc
pars <- spimalot::spim_fit_pars_load("parameters", region, "central",
                                     kernel_scaling)

pars <- simplify_transform(pars, "parameters", date)

## Fix all unused parameters (those not impacting fitting before the date parameter)
pars <- fix_unused_parameters(pars, date)

#Create the conversions functions for the parameters in order to fit in |R^d
#And attach them to the pars object
pars <- create_Rn2par(pars)

restart_date <- NULL

## This will probably want much more control, if we are to support
## using rrq etc to create a multinode job; some of that will depend a
## bit on the combination of multiregion and deterministic I think

#This sets up a lot of pmcmc controls, checks iterations are compatible etc.
control <- spimalot::spim_control(
  FALSE, chains, deterministic, date_restart = restart_date,
  n_mcmc = n_mcmc, burnin = burnin,
  compiled_compare = deterministic, adaptive_proposal = deterministic)


data_rtm <- read_csv("data/rtm.csv")
data_serology <- read_csv("data/serology.csv")

#This trims off dates and columns we don't want
data <- spim_data(
  date, region, model_type, data_rtm, data_serology,
  trim_deaths, trim_pillar2,
  full_data = FALSE)

#Some more particle filter setup, this object can then be "run" below
filter <- spimalot::spim_particle_filter(data, pars$mcmc,
                                         control$particle_filter,
                                         deterministic)

## To run the model at this point, we just need to run:
##
## > filter$run(pars$mcmc$model(pars$mcmc$initial()))

## This bit takes ages, of course
samples <- spimalot::spim_fit_run(pars, filter, control$pmcmc)

write.csv(samples$pars,file = 'samples.csv')

## This is the data set including series that we do not fit to, and
## with the full series of carehomes deaths.
data_full <- spim_data(
  date, region, model_type, data_rtm,
  data_serology, trim_deaths, trim_pillar2,
  full_data = TRUE)

## This is new, and used only in sorting out the final outputs. Some
## explanation would be useful.
data_inputs <- list(rtm = data_rtm,
                    full = data_full,
                    fitted = data)

dat <- spimalot::spim_fit_process(samples, pars, data_inputs,
                                  control$particle_filter)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat$fit, "outputs/fit.rds")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$fit$samples),
  width = 16, height = 9)
