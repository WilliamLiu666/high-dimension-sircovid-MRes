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

epsilon <- 0.001
L <- 1
N <- 5001
HMC_samples <- matrix(0,N,length(theta))
HMC_samples[1,] <- theta
M <- diag(1,length(theta))
invM <- M

## Run HMC with 5 blocks
for (j in 1:1){

  ## HMC for N iterations
  for (i in ((j-1)*N+2):(j*N+1)){
    print(i)
    HMC_samples[i,] <- HMC_parallel(RnPosterior, gradient_LP, epsilon, L, HMC_samples[i-1,], filter,filter2, pars, M, invM)
  }

  #Update variance matrix
  invM <- cov(HMC_samples[((j-1)*N+2):(j*N+1),])
  M <- solve(invM)
}


## save the output
write.csv(HMC_samples,file = sprintf('samples_%s_%s_acc4.csv',epsilon,L))


covmat <- read_csv("covmat.csv")
M <- as.matrix(covmat[,2:29])
invM <- diag(diag(M),length((theta)))
M <- solve(invM)

## HMC parameters
L <- 1
N <- 1001
ESS <- matrix(0,nrow = 6, ncol = 28)
acc.list <- rep(0,6)
epsilon <- 0.01

for (j in 1:2){
  entry <- (j*.01+2.89)*1e-7
  invM[7,7] <- entry
  M <- solve(invM)
  HMC_samples <- matrix(0,N,length(theta))
  HMC_samples[1,] <- theta
  for (i in 2:N){
    print(sprintf('%s th iter with 7th diagonal entry = %s .csv',i,entry))
    HMC_samples[i,] <- HMC_parallel(RnPosterior, gradient_LP, epsilon, L, HMC_samples[i-1,], filter,filter2, pars, M, invM)
  }
  write.csv(HMC_samples,file = sprintf('samples_7th_entry_%s.csv',entry))
  ESS[j,] <- effectiveSize(HMC_samples)
  acc.list[j] <- acc_rate(HMC_samples)
}

# for (j in 1:11){
#   entry <- (j*.1+2.4)*1e-7
#   HMC_samples <- read_csv(sprintf('samples_7th_entry_%s.csv',entry))
#   HMC_samples <- as.matrix(HMC_samples[,2:29])
#   ESS[j,] <- effectiveSize(HMC_samples)
#   acc.list[j] <- acc_rate(HMC_samples)
# }