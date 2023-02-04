library(MASS)
library(coda)
simulate_calculate_vaccination_new <- function(state, pars, region) {
  n_groups <- sircovid:::lancelot_n_groups()
  
  ## output the cumulative transitions between vaccine strata
  ## by age / vaccine stratum / region / over time
  n_vaccinated <- apply(state[startsWith(rownames(state),"cu"), , , , drop = FALSE],
                        c(1, 3, 4), mean)
  n_strata <- nrow(n_vaccinated) / n_groups
  n_vaccinated <-
    mcstate::array_reshape(n_vaccinated, 1L, c(n_groups, n_strata))
  
  
  # Output number of first, second and booster doses
  
  idx_doses <- c("first_dose" = 1, "second_dose" = 2, "waned_dose" = 3, "booster_dose" = 4)
  doses <- n_vaccinated[, idx_doses, , ,drop = FALSE]
  dimnames(doses)[2:3] <- list(names(idx_doses), region)
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")
  
   abind_quiet(doses, doses_inc, along = 2)
  
}

#functions to pass parameter from Rn to par and vice-versa
#This builds the function and use eval-parse to create it
create_Rn2par <- function(pars){
  
  #use this list as the reference for the list of names
  #the Rn vector is matching the order of this list
  par_names <- names(pars$mcmc$initial())
  
  i <- match(par_names, pars$info$name)
  pars_min <- pars$info$min[i]
  pars_max <- pars$info$max[i]
  names(pars_min) <- names(pars_max) <- par_names
  
  pars$par2Rn <- function(x) log((x - pars_min) / (pars_max - x)) 
  pars$Rn2par <- function(x) (pars_min + pars_max * exp(x)) / (1 + exp(x))
  
  pars
}

#This function takes a R^n function and calculate the posterior of our model
RnPosterior <- function(theta, filter, pars){
  p <- pars$Rn2par(theta)
  LL_theta <- filter$run(pars = pars$mcmc$model(p))
  LP_theta <- pars$mcmc$prior(p)
  return(LL_theta+LP_theta)
}

#This function evaluate the posterior at multiple point in the parameter space
calculate_posterior_map <- function(parameter_samples, filter, pars){
  apply(parameter_samples, 2, function(x)
    RnPosterior(x, filter, pars))
}

#Calculate the gradient
gradient_LP <- function(theta, filter, pars, eps = 1e-4){
  n <- length(theta)
  theta_h <- diag(eps,n) + matrix(rep(theta,n),ncol = n)
  
  LP_theta <- RnPosterior(theta, filter, pars)
  LP_h <- calculate_posterior_map(theta_h, filter, pars)
  
  list(LP = LP_theta,
       grad_LP = (LP_h-LP_theta)/eps)
}

gradient_LP_parallel <- function(theta, filter, pars, eps = 1e-4){
  n <- length(theta)
  ## first column will be theta, the following n columns will correspond to
  ## one parameter perturbed
  theta_map <- cbind(rep(0, n), diag(eps, n)) + matrix(rep(theta, n + 1), ncol = n + 1)

  ## calculate posterior across all columns in theta_map
  LP <- calculate_posterior_map_parallel(theta_map, filter, pars)
  
  ## first value corresponds to theta
  LP_theta <- LP[1]
  ## the rest used to calculate gradient estimate
  LP_h <- LP[-1]
  names(LP_h) <- names(theta)
  
  # #one parameter perturbed
  # theta_map <- cbind(rep(0, n), diag(-eps, n)) + matrix(rep(theta, n + 1), ncol = n + 1)
  #
  # #calculate posterior across all columns in theta_map
  # LP <- calculate_posterior_map_parallel(theta_map, filter, pars)
  #
  # first value corresponds to theta
  # LP_theta <- LP[1]
  # #the rest used to calculate gradient estimate
  # LP_h2 <- LP[-1]
  # names(LP_h2) <- names(theta)
  
  list(LP = LP_theta,
       grad_LP = (LP_h - LP_theta) / (eps))
}

#This function evaluate the posterior at multiple point in the parameter space
calculate_posterior_map_parallel <- function(parameter_samples, filter, pars){
  ## transform from Rn to parameter space
  p <- apply(parameter_samples, 2, pars$Rn2par)
  
  ## transform to odin parameters
  model_pars <- apply(p, 2, pars$mcmc$model)
  
  ## calculate log-likelihoods in parallel
  LL_theta <- filter$run(model_pars)
  ## calculate log-priors
  LPr_theta <- apply(p, 2, pars$mcmc$prior)
  
  LL_theta + LPr_theta
  
}

resize_filter <- function(filter, n_parameters, n_threads) {
  inputs <- filter$inputs()
  inputs$n_parameters <- n_parameters
  inputs$n_threads <- n_threads
  ## Note that this is not yet exported from mcstate, I'll look into
  ## that for you later (mrc-3868 for my reference)
  mcstate:::particle_filter_from_inputs(inputs)
}


## This function will fix any unused parameters that do not impact fitting 
## before the end date
fix_unused_parameters <- function(pars, date) {
  
  ## Automatically detect which betas to fix
  ## We need to keep all betas up to the first one with date greater than or 
  ## equal to the date parameter
  i <- max(which(pars$base$beta_date < sircovid::sircovid_date(date)))
  beta_fixed <- setdiff(pars$base$beta_names, sprintf("beta%d", seq_len(i + 1)))
  
  ## Now we will fix other parameters that have no impact before the date
  ## parameter
  
  ## Note firstly that the following parameters are required whatever
  ## the date parameter is:
  ## "alpha_D", "alpha_death_hosp", "alpha_H", "eps", "m_CHR", "m_CHW", 
  ## "p_G_D", "p_G_D_CHR", "p_H", "p_H_D", "p_ICU", "p_ICU_D", "p_W_D",
  ## "start_date"
  
  ## Now we declare the date from which these parameters have an impact
  pars_dates <- list(
    ## Various changepoint parameters. We must include them the day after the
    ## previous changepoint
    mu_D = "2020-04-02",
    mu_D_2 = "2020-10-02",
    mu_gamma_H = "2020-12-02",
    mu_gamma_H_2 = "2021-01-02", 
    mu_gamma_H_3 = "2021-03-02", 
    mu_gamma_H_4 = "2021-06-02",
    p_H_2 = "2020-10-02",
    p_H_3 = "2020-12-16", 
    p_ICU_2 = "2020-04-02",
    
    ## Pillar 2 parameters - we start fitting pillar 2 from 2020-06-18
    ## This is a Thursday, so first weekend day is 2020-06-20
    p_NC = "2020-06-18",
    p_NC_weekend = "2020-06-20", 
    rho_pillar2_tests = "2020-06-18",
    
    ## delta parameters
    rel_p_D_delta = "2021-03-08",
    rel_p_H_delta = "2021-03-08",
    rel_p_ICU_delta = "2021-03-08",
    seed_date_delta = "2021-03-08",
    ta_delta = "2021-03-08"
  )
  
  ## Fix parameters
  fixed <- c(beta_fixed, names(pars_dates)[pars_dates > date])
  pars$mcmc <- pars$mcmc$fix(pars$mcmc$initial()[fixed])
  
  pars
}


HMC_parallel <- function (RnPosterior, gradient_LP, epsilon, L, current_q, filter,filter2, pars, M, invM)
{
  q = current_q
  
  p = mvrnorm(1,rep(0,length(q)),M) # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  p = p - epsilon * -gradient_LP_parallel(q, filter2, pars)$grad_LP / 2
  
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * invM%*%p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * -gradient_LP_parallel(q, filter2, pars)$grad_LP
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * -gradient_LP_parallel(q, filter2, pars)$grad_LP / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = -RnPosterior(current_q, filter, pars)
  current_K = current_p %*% invM %*% current_p /2
  proposed_U = -RnPosterior(q, filter, pars)
  proposed_K = p %*% invM %*% p /2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
  
}

HMC <- function (RnPosterior, gradient_LP, epsilon, L, current_q, filter, pars)
{
  q = current_q
  
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  p = p - epsilon * -gradient_LP(q, filter, pars)$grad_LP / 2
  
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * -gradient_LP(q, filter, pars)$grad_LP
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * -gradient_LP(q, filter, pars)$grad_LP / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = -RnPosterior(current_q, filter, pars)
  current_K = sum(current_p^2) / 2
  proposed_U = -RnPosterior(q, filter, pars)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
  
}

acc_rate <- function(x){
  num <- dim(x)[1]-1
  count <- 0
  for (i in 1:num){
    if (x[i+1]==x[i]){
      count<-count+1
    }
  }
  return((num-count)/num)
}


simplify_transform <- function(pars, path, date) {
  
  e <- new.env()
  sys.source(file.path(path, "transform.R"), e)
  
  make_transform <- e$make_transform
  pars$transform <- make_transform(pars$base, date)
  
  pars$mcmc <- spimalot:::spim_pars_mcmc_single(pars$info, pars$prior, 
                                                pars$proposal, pars$transform)
  
  pars$base$epoch_dates <-
    pars$base$epoch_dates[pars$base$epoch_dates <= sircovid_date(date)]

  pars
}

