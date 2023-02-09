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


filter$run(pars$mcmc$model(pars$mcmc$initial()))
n_threads <- spimalot::spim_control_cores()
n_pars <- length(pars$mcmc$initial())
## We use n_pars + 1 here as we calculate the log-likelihood at theta and
## then also at a perturbation in each dimension

filter2 <- resize_filter(filter, n_pars + 1, n_threads)

## load parameters
theta <- read_csv("theta.csv")
theta <- as.array(theta$x)
invM <- read.csv("invM.csv")
invM <- as.matrix(invM[,2:29])
M <- solve(invM)


HMC_samples <- matrix(0,N+1,length(theta))
HMC_samples[1,] <- theta

if (compare == TRUE){
  acc.list <- matrix(0,nrow =N, ncol = 6)
}else{
  acc.list <- rep(0,N)
}


## HMC for N iterations
for (i in 1:N){
  ind <- i+1
  if (ind%%10 == 0){
    print(ind)
  }
  result <- HMC_parallel(RnPosterior, gradient_LP_parallel, epsilon, L, HMC_samples[ind-1,], filter,filter2, pars, M, invM, compare = compare)
  HMC_samples[ind,] <- result$q
  if (compare == TRUE){
    acc.list[ind-1,] <- result$acc.list
  }
  else{
    acc.list[ind-1] <- result$acc.list
  }
}



message("Saving results")
dir.create("outputs", FALSE, TRUE)


message("Saving plots")

