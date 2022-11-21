
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
  
  #Start of the function
  Rn2par_text <- "Rn2par <- function(x){"
  par2Rn_text <- "par2Rn <- function(x){"
  
  #Create a name 0 valued Named num vector in the par space
  #and not named in the Rn space
  Rn2par_text <- paste(Rn2par_text, "res <- rep(0,",length(par_names), ")\n", sep="")
  Rn2par_text <- paste0(Rn2par_text, "names(res) <- c(\"",
                        paste(par_names, collapse="\",\""),"\")\n")
  par2Rn_text <- paste(par2Rn_text, "res <- rep(0,",length(par_names), ")\n", sep="")
  
  #Loop over the parameter names
  for(i in seq_along(par_names)){
    
    #detect if the parameter is from a beta, gamma or undefined distribution
    dist_form <- pars$prior[pars$prior$name==par_names[i],]$type
    
    if(dist_form == "gamma"){
      Rn2par_text <- paste0(Rn2par_text, "res[\"",
                            par_names[i], "\"] <- exp(x[",i,"])\n")
      par2Rn_text <- paste0(par2Rn_text, "res[",i,"] <- log(x[\"",
                            par_names[i], "\"","])\n")
      
    }
    
    if(dist_form == "beta"){
      Rn2par_text <- paste0(Rn2par_text, "res[\"",
                            par_names[i], "\"] <- exp(x[",i,"])/(exp(x[",i,"])+1)\n")
      par2Rn_text <- paste0(par2Rn_text, "res[",i,"] <- log(x[\"",
                            par_names[i], "\"","]/(1-x[\"",
                            par_names[i], "\"","]))\n")
    }
    
    if(dist_form == "null")
    {
      Rn2par_text <- paste0(Rn2par_text, "res[\"",
                            par_names[i], "\"] <- x[",i,"]\n")
      par2Rn_text <- paste0(par2Rn_text, "res[",i,"] <- x[\"",
                            par_names[i], "\"]\n")
    }
  }
  #Return the whole vector
  Rn2par_text <- paste0(Rn2par_text, "res \n}")
  par2Rn_text <- paste0(par2Rn_text, "res \n}")
  return(list(Rn2par = Rn2par_text, par2Rn = par2Rn_text))
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
gradient_LP <- function(theta, pars, filter, eps = 1e-4){
  n <- length(theta)
  theta_h <- diag(eps,n) + matrix(rep(theta,n),ncol = n)
  
  LP_theta <- RnPosterior(theta, filter, pars)
  LP_h <- calculate_posterior_map(theta_h, filter, pars)
  
  list(LP = LP_theta,
       grad_LP = (LP_h-LP_theta)/eps)
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
