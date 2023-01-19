source("global_util.R")

version_check("sircovid", "0.14.9")
version_check("spimalot", "0.8.16")

## Can change the end date for the fitting period here. These tasks were setup
## with an end date of 2021-09-13, but an earlier date can be used
date <- "2020-09-13"
model_type <- "BB"

#Some data streams are unreliable for the last few days. Define how many days here.
trim_deaths <- 4
trim_pillar2 <- 5

## MCMC control (only applies if short_run = FALSE)
burnin <- 8000
n_mcmc <- 10000
chains <- 1
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
# grad <- gradient_LP_parallel(theta, filter2, pars)

## Using MCMC to find a better initial value.
samples <- spimalot::spim_fit_run(pars, filter, control$pmcmc)
theta <- samples$pars[1000,]
theta <- as.vector(theta)
theta <- pars$par2Rn(theta)

## HMC parameters
epsilon <- 0.0015
L <- 1
N <- 5001
HMC_samples <- matrix(0,N,length(theta))
HMC_samples[1,] <- theta
library(MASS)
M <- diag(1,length(theta))
invM <- M

## Run HMC with 5 blocks
for (j in 1:5){
  
  ## HMC for N iterations
  for (i in ((j-1)*N+2):(j*N+1)){
    print(i)
    HMC_samples[i,] <- HMC_parallel(RnPosterior, gradient_LP, epsilon, L, HMC_samples[i-1,], filter,filter2, pars, M, invM)
  }
  
  #Update variance matrix
  M <- cov(HMC_samples[((j-1)*N+2):(j*N+1),])
  invM <- diag(diag(M),length((theta)))
  M <- solve(invM)
}


## save the output
write.csv(HMC_samples,file = sprintf('samples_%s_%s.csv',epsilon,L))

pdf(sprintf('HMC_epsilon_%s_L_%s_noburnin.pdf',epsilon,L))
par(mfrow=c(3,4))
for (i in 1:6){
  acf(HMC_samples[,i],main = sprintf('the %s th feature',i))
  plot(HMC_samples[,i],main = sprintf('the %s th feature',i))
}
for (i in 7:12){
  acf(HMC_samples[,i],main = sprintf('the %s th feature',i))
  plot(HMC_samples[,i],main = sprintf('the %s th feature',i))
}
for (i in 13:18){
  acf(HMC_samples[,i],main = sprintf('the %s th feature',i))
  plot(HMC_samples[,i],main = sprintf('the %s th feature',i))
}
for (i in 19:24){
  acf(HMC_samples[,i],main = sprintf('the %s th feature',i))
  plot(HMC_samples[,i],main = sprintf('the %s th feature',i))
}
for (i in 25:28){
  acf(HMC_samples[,i],main = sprintf('the %s th feature',i))
  plot(HMC_samples[,i],main = sprintf('the %s th feature',i))
}
dev.off()



















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
dat$fit$simulate$n_doses <-
simulate_calculate_vaccination_new(dat$fit$simulate$state, pars, region)

dir.create("outputs", FALSE, TRUE)
saveRDS(dat$fit, "outputs/fit.rds")
saveRDS(dat$restart, "outputs/restart.rds")

message("Creating plots")
write_pdf(
  "outputs/pmcmc_traceplots.pdf",
  spimalot::spim_plot_fit_traces(dat$fit$samples),
  width = 16, height = 9)
