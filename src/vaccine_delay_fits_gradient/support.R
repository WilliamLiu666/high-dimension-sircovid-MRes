
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
RnPosterior <- function(theta){
  p <- Rn2par(theta)
  LL_theta <- filter$run(pars = pars$mcmc$model(p))
  LP_theta <- pars$mcmc$prior(p)
  return(LL_theta+LP_theta)
}

#This function evaluate the posterior at multiple point in the parameter space
calculate_posterior_map <- function(parameter_samples){
  apply(parameter_samples, 2, function(x)
    RnPosterior(x))
}

#Calculate the gradient
gradient_LP <- function(theta, eps = 1e-4){
  n <- length(theta)
  theta_h <- diag(eps,n) + matrix(rep(theta,n),ncol = n)
  
  LP_theta <- RnPosterior(theta)
  LP_h <- calculate_posterior_map(theta_h)
  
  list(LP = LP_theta,
       grad_LP = (LP_h-LP_theta)/eps)
}