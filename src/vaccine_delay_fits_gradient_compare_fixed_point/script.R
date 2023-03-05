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

grad <- array(0, c(12,6,28))
axis.x <- rep(0,12)
for (i in 2:13){
  eps = 10^-i
  axis.x[i-1] <- eps
  grad[i-1,,] <- gradient_LP_parallel(theta,filter2,pars,method = 'all', eps = eps)
}


# true <- 0
# dimension <- 2
# df <- data.frame(central_1st_order = log(abs(grad[2:11,1,dimension]-true)),central_2nd_order = log(abs(grad[2:11,2,dimension]-true)),
#                  forward_1st_order = log(abs(grad[2:11,3,dimension]-true)),forward_2nd_order = log(abs(grad[2:11,4,dimension]-true)),
#                  backward_1st_order = log(abs(grad[2:11,5,dimension]-true)),backward_2nd_order = log(abs(grad[2:11,6,dimension]-true)),
#                  axis.x = axis.x[2:11])

df <- data.frame(central_1st_order = (rowSums(grad[2:11,1,]^2))^(1/2),central_2nd_order = (rowSums(grad[2:11,2,]^2))^(1/2),
                 forward_1st_order = (rowSums(grad[2:11,3,]^2))^(1/2),forward_2nd_order = (rowSums(grad[2:11,4,]^2))^(1/2),
                 backward_1st_order = (rowSums(grad[2:11,5,]^2))^(1/2),backward_2nd_order = (rowSums(grad[2:11,6,]^2))^(1/2),
                 axis.x = axis.x[2:11])


# message("Saving results")
# dir.create("outputs", FALSE, TRUE)
# write.csv(grad,'outputs/gradient.csv')


# message("Saving plots")
library(ggplot2)
p <- ggplot(df,aes(x = axis.x))+
  geom_line(aes(y = central_1st_order,color = "central 1st order"),lwd=1.2)+
  geom_line(aes(y = central_2nd_order,color = "central 2nd order"),lwd=1.2)+
  geom_line(aes(y = forward_1st_order,color = "forward 1st order"),lwd=1.2)+
  geom_line(aes(y = forward_2nd_order,color = "forward 2nd order"),lwd=1.2)+
  geom_line(aes(y = backward_1st_order,color = "backward 1st order"),lwd=1.2)+
  geom_line(aes(y = backward_2nd_order,color = "backward 2nd order"),lwd=1.2)+
  scale_x_continuous(trans = 'log10')+
  labs(x = "h" , y =  "2-norm of gradient",
                                           color = "Legend")

# pdf('outputs/gradient_plot.pdf')
p
dev.off()

