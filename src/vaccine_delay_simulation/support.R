spim_simulate_local <- function(pars, combined, rt_critical_dates, date_end,
                                output, n_threads, seed) {
  ## This ensures a reproducible seed (starting from R's RNG seed)
  ## regardless of how many threads are used, and if we distribute
  ## over multiple nodes later.
  seed <- dust::dust_rng_distributed_state(seed, n_nodes = length(pars))
  lapply(seq_along(pars), spim_simulate, pars, combined, rt_critical_dates,
         date_end, output, n_threads, seed)
}


spim_simulate <- function(i, pars, combined, rt_critical_dates, date_end,
                          output, n_threads, seed) {
  #browser()
  el <- pars[[i]]
  message(sprintf("-----\nRunning scenario %d / %d", i, length(pars)))
  time <- system.time(
    ret <- spim_simulate_one(el, combined, rt_critical_dates, date_end, output,
                             n_threads, seed[[i]]))
  message(sprintf("Finished scenario %d in %2.1f s", i, time[["elapsed"]]))
  ret
}


spim_simulate_one <- function(pars, combined, rt_critical_dates, date_end, output,
                              n_threads, seed) {

  regions <- names(pars)

  step_start <- combined$step
  date_start <- sircovid::sircovid_date(combined$date)
  date_end <- sircovid::sircovid_date(date_end)
  dates <- seq(date_start, date_end)
  steps <- dates * combined$steps_per_day
  step_end <- last(steps)
  state_start <- combined$state

  pars <- unlist(pars, FALSE, FALSE)
  attributes(pars) <- attributes(combined$pars) # dim and dimnames

  message("Creating dust object")
  obj <- sircovid::lancelot$new(pars, step_start, NULL, pars_multi = TRUE,
                                n_threads = n_threads, seed = seed)

  #----
  info <- obj$info()[[1]]
  index <- spimalot:::simulate_index(info,
                                     output$keep,
                                     output$vaccination,
                                     multistrain = TRUE)

  ## TODO: this is a hack for now, we'll need to improve this in
  ## spimalot, but that depends a bit on what we do with this.
  cum_infections_per_strain <- info$index$cum_infections_per_strain

  names(cum_infections_per_strain) <-
    sprintf("cum_infections_per_strain_%s",
            seq_along(cum_infections_per_strain))
  index$run <- c(index$run, cum_infections_per_strain)
  index$cum_infections_per_strain <- cum_infections_per_strain

  ## Zero infections before we start, so that counts are since the
  ## beginning of the simulation.  Probably we need to add the seed
  ## here to be completely accurate though.
  state_start[cum_infections_per_strain, , ] <- 0

  obj$update_state(state = state_start)
  obj$set_index(index$run)
  message("Simulating!")
  state <- obj$simulate(steps)
  dimnames(state)[[3]] <- regions

  ## So strain 1 is Delta
  ##    strain 2 is Omicron
  ##    strain 3 is Omicron after Delta
  ##    strain 4 is Delta after Omicron
  ## Take the sum of strains 2:3 over all to get the fraction omicron
  idx <- names(index$cum_infections_per_strain)
  n_strain_1 <- apply(state[idx[c(1, 4)], , , ], 2:4, sum)
  n_strain_2 <- apply(state[idx[c(2, 3)], , , ], 2:4, sum)

  message("Adding summary statistics")
  ret <- list(
    date = dates,
    summary_state = spimalot:::create_summary_state(
      state, output$keep, dates))

  if (output$time_series) {
    ## TODO: this is a hack to create some space to put the fractions
    ## into the output. We could easily expand this to collect all 4
    ## metastrains.
    ret$state <- state[c(output$keep, output$keep[c(1, 1)]), , , ]
    nr <- nrow(ret$state)
    ret$state[nr - 1, , , ] <- n_strain_1
    ret$state[nr    , , , ] <- n_strain_2
    rownames(ret$state)[c(nr - 1, nr)] <- c("n_strain_1", "n_strain_2")
  }

  if (output$rt) {
    message("Calculating Rt")
    rt <- spimalot:::simulate_rt(
      steps,
      state[names(index$S), , , ],
      pars,
      rt_critical_dates,
      TRUE,
      state[names(index$R), , , ],
      state[names(index$prob_strain), , , ],
      no_seeding = TRUE,
      prop_voc = NULL,
      weight_Rt = FALSE,
      rt_type = output$rt_type
    )

    if (output$rt_weighted) {
      rt_weighted <- spimalot:::simulate_rt(
        steps,
        state[names(index$S), , , ],
        pars,
        rt_critical_dates,
        TRUE,
        state[names(index$R), , , ],
        state[names(index$prob_strain), , , ],
        no_seeding = TRUE,
        prop_voc = NULL,
        weight_Rt = TRUE,
        rt_type = output$rt_type)

      rt_join <- function(rt, weighted_rt) {
        x <- abind_quiet(rt, weighted_rt, along = 4)
        dimnames(x)[[4]] <- c("strain_1", "strain_2", "both")
        x
      }

      # combined weighted and strain specific outputs
      rt <- Map(rt_join, rt, rt_weighted)
    }
    ret <- c(ret, rt)
  }

  if (output$state_by_age) {
    ret$state_by_age <- spimalot:::simulate_extract_age_class_state(
      state, index)
    ## TODO: this was previously set within spimalot by looking at the
    ## number of vaccination categories, but that is unsafe!  Adding
    ## some here, but this will need changing somewhere else.
    for (i in seq_along(ret$state_by_age)) {
      colnames(ret$state_by_age[[i]]) <-
        c("unvaccinated", "partial_protection", "full_protection",
          "waned_protection")
      #, "booster") ##Removed because we're not using boosters
    }
  }
  
  n_strain <- 4

  if (output$vaccination) {
    message("Computing vaccination")
    ## TODO: it will be necessary to update the functions
    ## calculate_n_protected/simulate_calculate_vaccination so that
    ## they naturally work in terms of sircovid parameters as that is
    ## what we have at this point.

    ## If we take the first parameters element (pars[[1]]) that has
    ## all the parameters for computing vaccine efficacy and will know
    ## things like if boosters are involved, how many strains there
    ## are etc.  This needs chasing through a lot of functions in
    ## spimalot and someone from the epi side should make sure it's
    ## done correctly.
    ##
    ## NOTE: assumption here is that vaccine efficacy is constant over
    ## a block of parameters (i.e., not varying by region or by
    ## particle), which is a reasonable assumption at the moment.
    vaccination_summary <- simulate_calculate_vaccination(
      state, index, pars[[1]])
    ret <- c(ret, vaccination_summary)
  }

  ret
}


prepare_parameter_update <- function(dat, combined, control) {
  regions <- control$regions
  date_end <- sircovid::sircovid_date(control$date_end)
  lapply(dat, function(x)
    set_names(
      lapply(regions, prepare_parameter_update1, x, combined, date_end),
      regions))
}


## TODO: move to utils
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


prepare_parameter_update1 <- function(region, pars_simulation, combined,
                                      date_end) {
  #browser()
  pars <- combined$pars[, region]
  pars1 <- pars[[1]]
  ## TODO: It's not clear what will or will not work here as we change
  ## the number of strains around really.  Ideally we'll do all the
  ## munging elsewhere eventually.
  stopifnot(pars1$n_strains == 4)

  ## TODO: think we want to move away from doing this and rather extract all
  ## baseline parameters we can possibly need from the fits in
  ## rtm_parameters_simulate
  baseline <- combined$base[[region]]
##TODO: REMOVE BOOSTER STUFF
  booster_eligibility <-
    sircovid::vaccine_eligibility(pars_simulation$booster_eligibility_min_age)

  ## We can create a list of things to update this way:
  vaccine_schedule <- compute_vaccine_schedule(
    region, baseline, pars_simulation, date_end, pars1,
    booster_eligibility)

  vaccine_schedule_effect <-
    shift_doses(vaccine_schedule, baseline$vaccine_days_to_effect)

  #THOM EDIT:

  #Here, if the orderly parameter isn't FALSE, we set a brand new rel_severity_alpha_delta
  if(!isFALSE(vaccine_efficacy_alteration)){
    #Load data
    vaccine_efficacy_alpha <- read_csv("vaccine_data/vaccine_efficacy_alpha.csv")
    vaccine_efficacy_delta <- read_csv("vaccine_data/vaccine_efficacy_delta.csv")
    data_vaccination <- read_csv("vaccine_data/data_vaccination.csv")
    
    #An option if we want the very targetted changes:
    if(vaccine_efficacy_alteration == 'targetted'){
      vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose != 1 & vaccine_efficacy_alpha$type %in% c("disease", "infection")] <- vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose != 1 & vaccine_efficacy_alpha$type %in% c("disease", "infection")] - 0.1
      vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose != 1 & vaccine_efficacy_delta$type %in% c("disease", "infection")] <- vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose != 1 & vaccine_efficacy_delta$type %in% c("disease", "infection")] - 0.1
      
    }else if(vaccine_efficacy_alteration == "first_dose"){
    #This is the sensitivity where all first dose VEs are increase by 0.1
      vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 1] <- vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 1] + 0.1
      vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 1] <- vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 1] + 0.1
      
      #Don't let 1dose VE be greater than 2dose VE. Cap any increases at 2 dose
      vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 1] <- pmin(vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 2],
                                                                              vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 1])
      
      #Don't let 1dose VE be greater than 2dose VE. Cap any increases at 2 dose
      vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 1] <- pmin(vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 2],
                                                                              vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 1])
      
    }else{ 
    #Now make the alteration to values
      vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose != 1] <- vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose != 1] + vaccine_efficacy_alteration
      vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose != 1] <- vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose != 1] + vaccine_efficacy_alteration
      
      #Don't let 1dose VE be greater than 2dose VE. Cap any decreases at 1 dose
      vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 2] <- pmax(vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 2],
                                                                              vaccine_efficacy_delta$central[vaccine_efficacy_delta$dose == 1])
      
      vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 2] <- pmax(vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 2],
                                                                               vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$dose == 1])
      
    }
    #Don't let values be less than 0 or above 1
    vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$central < 0] <- 0
    vaccine_efficacy_delta$central[vaccine_efficacy_delta$central < 0] <- 0
    vaccine_efficacy_alpha$central[vaccine_efficacy_alpha$central > 1] <- 1
    vaccine_efficacy_delta$central[vaccine_efficacy_delta$central > 1] <- 1
    

    
    # Data on vaccine doses by age, date and type of vaccine
    vacc_doses_by_age_date_vaccine <-
      data_vaccination %>%
      dplyr::mutate(age_band_min = replace(age_band_min, age_band_min == 16, 15)) %>%
      dplyr::filter(!is.na(age_band_min)) %>%
      dplyr::group_by(vaccine, age_band_min) %>%
      dplyr::summarise(n = sum(first_dose, na.rm = TRUE)) %>%
      dplyr::group_by(age_band_min) %>%
      dplyr::mutate(freq = n / sum(n, na.rm = TRUE)) %>%
      tidyr::pivot_wider(id_cols = c(age_band_min), values_from = freq, names_from = vaccine)
    
    ## check all proportions sum to 1
    stopifnot(
      all(
        abs(rowSums(vacc_doses_by_age_date_vaccine[3:17, -1], na.rm = TRUE) - 1) < 1e-8
      )
    )
    
    prop_pfizer <- vacc_doses_by_age_date_vaccine$Pfizer + vacc_doses_by_age_date_vaccine$Moderna
    
    ## set %PF in CHW / CHR <-  80+
    prop_pfizer <- c(prop_pfizer, rep(prop_pfizer[17], 2))
    # 100% PF in 0-10
    prop_pfizer[sircovid:::sircovid_age_bins()$end < 10] <- 1
    #browser()
    ## Note that we do not change efficacy to allow for differences in future uptake
    average_vacc_efficacy_alpha <-
      calculate_average_vacc_efficacy(vaccine_efficacy_alpha, prop_pfizer)
    
    ## format that is readable by sircovid
    ##
    ## TODO: Simplify the csv file too and the loading here as we will
    ## only use a central estimate now (these apply only to alpha)
    rel_severity_alpha <- lapply(average_vacc_efficacy_alpha, function(e)
      get_vaccine_conditional_prob(e$death,
                                   e$severe_disease,
                                   e$disease,
                                   e$infection,
                                   e$transmission))$central
    
    
    ## This all can be refactored as it's not actually doing anything
    ## nearly as complicated here as we are doing as this was
    ## previously much more general.
    average_vacc_efficacy_delta <-
      calculate_average_vacc_efficacy(vaccine_efficacy_delta, prop_pfizer)
    rel_severity_delta <- lapply(
      average_vacc_efficacy_delta, function(e)
        get_vaccine_conditional_prob(e$death, e$severe_disease, e$disease,
                                     e$infection, e$transmission))$central
    
    rel_severity_alpha_delta <-
      rel_severity_strains(rel_severity_alpha,
                           rel_severity_delta)
    
    pars_simulation$strain_vaccine_efficacy <- rel_severity_alpha_delta
    
  }
  
  ## NOTE: this is not actually allowed above, but I think it would be
  ## easy enough here.
  vaccine_efficacy <-
    pars_simulation$strain_vaccine_efficacy %||%
    baseline$rel_severity_alpha_delta
  #Remove all fifth columns of 3rd dim
  for(i in seq_along(vaccine_efficacy)) {
    
    vaccine_efficacy[[i]] <- vaccine_efficacy[[i]][, , -5]
    
  }
  
  ##Here we can overwrite the mean time to waning
  
  if(!isFALSE(mean_waning_time)){
    pars_simulation$vaccine_waning_days <- mean_waning_time
  }
  
  
  ## Update vaccine progression to use the latest time to VE waning
  vaccine_progression_rate <- baseline$vaccine_progression_rate
  if (!is.null(pars_simulation$vaccine_waning_days)) {
    vaccine_progression_rate[baseline$vaccine_index_dose2 + 1] <-
      1 / pars_simulation$vaccine_waning_days
  }
  #I am now removing the fifth and final entry of vaccine_progression_rate
  vaccine_progression_rate <- vaccine_progression_rate[-5]

  ## TODO: If we add boosters by inflation, this is going to be
  ## problematic. In particular we should make sure that we update the
  ## baseline in the inflated group.
  ##
  ## TODO: would be nicer to have a new function like
  ## sircovid::update_pars_vaccination_schedule that checks that the
  ## new vaccination schedule is compatible with the rest of the
  ## parameters rather than building it up from scratch.
  ##
  ## TODO: There is a bug here! There is post-processing that happens
  ## in sircovid to update the protection against death and that needs
  ## to be done here.
  ## https://github.com/mrc-ide/sircovid/blob/4fd526/R/lancelot.R#L509-L524
  update <- sircovid:::lancelot_parameters_vaccination(
    pars1$N_tot,
    pars1$dt,
    rel_susceptibility = vaccine_efficacy$rel_susceptibility,
    rel_p_sympt = vaccine_efficacy$rel_p_sympt,
    rel_p_hosp_if_sympt = vaccine_efficacy$rel_p_hosp_if_sympt,
    rel_p_death = vaccine_efficacy$rel_p_death,
    rel_infectivity = vaccine_efficacy$rel_infectivity,
    vaccine_schedule = vaccine_schedule_effect,
    vaccine_index_dose2 = baseline$vaccine_index_dose2,
    #vaccine_index_booster = baseline$vaccine_index_booster,
    vaccine_index_booster = NULL,
    vaccine_progression_rate = vaccine_progression_rate,
    n_strains = pars1$n_strains,
    n_doses = vaccine_schedule$n_doses)

  ## This is the only one that will need special treatment later as we
  ## modify each trajectory separately.
  beta_step <- pars_simulation$beta_step[, region, ]

  ## TODO: we will probably need to alter efficacy here from the
  ## baseline, particularly for boosters.
  ##
  ## TODO: This is a hack, previously we went through the (also
  ## unexported) function sircovid:::lancelot_parameters_strain.  This
  ## converts between the two-strain representation of parameters to 4
  ## and that will likely change as we move to a three-strain version.

  update$cross_immunity <-
    pars_simulation$strain_cross_immunity %||% baseline$cross_immunity

  if (!is.null(pars_simulation$serial_interval_scale)) {
    update$rel_gamma_E <-
      1 / sircovid:::mirror_strain(pars_simulation$serial_interval_scale)
  }

  ## beta_step is different as we vary it for each trajectory
  for (i in seq_along(pars)) {
    pars[[i]][names(update)] <- update
    pars[[i]]$beta_step <- beta_step[i, ]
  }

  pars
}


calculate_final_rt <- function(combined, parameters) {
  # Restart Rt used to prepare future Rt
  # England's Rt is used for every region in part of the process
  # but the ratio is used for each region independently
  ## TODO: should we do this pre-processing in the parameters task?
  england_rt <-
    apply(combined$simulate$Rt_general[, "england", , ], c(2, 3), mean)
  england_eff_rt <-
    apply(combined$simulate$eff_Rt_general[, "england", , ], c(2, 3), mean)

  ## Get last beta date from the combined object
  last_beta_date <- tail(combined$info$london$beta_date, 1)
  current_dates <- seq(last_beta_date, sircovid::sircovid_date(date))
  beta_mult_seasonality <- spimalot::spim_beta_mult_seasonality(
    current_dates, parameters$rt_seasonality_date_peak,
    parameters$rt_seasonality_amplitude
  )
  ## Based on this being *england* Rt, we only care if england's
  ## schools are open or closed. This is constant across regions, so
  ## pick one arbitrarily.
  beta_mult_schools <- spimalot::spim_beta_mult_schools(
    current_dates,
    control$parameters$rt_schools_schedule,
    control$parameters$rt_schools_modifier,
    "london")

  beta_mult <- beta_mult_seasonality * beta_mult_schools

  ## We manually remove the implicit seasonality from the fits when setting
  ##  Rt as this is then explicitly added when setting up beta step.
  avg_rt <- colMeans(england_rt[
    combined$simulate$date %in% current_dates, , drop = FALSE] / beta_mult)
  avg_eff_rt <- colMeans(england_eff_rt[
    combined$simulate$date %in% current_dates, , drop = FALSE] / beta_mult)

  data.frame(
    strain = c(names(avg_rt), names(avg_eff_rt)),
    rt_type = c(rep("rt", 3), rep("rt_eff", 3)),
    rt_value = c(avg_rt, avg_eff_rt),
    stringsAsFactors = FALSE
  )
}

## As strain transmission doesn't seem used it makes sense to abstract this
##  calculation to a new function that can be called only when required
## Function returns both the transmission advantage from the fits, and either:
##  1) strain transmission advantage from parameters;
##  2) c(1, advantage from fits)
calculate_strain_transmission_advantage <- function(rt, parameters) {
  fitted_s2s1_advantage <-
    rt[rt$strain == "strain_2" & rt$rt_type == "rt", "value"] /
      rt[rt$strain == "strain_1" & rt$rt_type == "rt", "value"]

  list(
    fitted_strain2_strain1_transmission_advantage = fitted_s2s1_advantage,
    strain_transmission = parameters$strain_transmission %||%
      c(1, fitted_s2s1_advantage)
  )
}

compute_vaccine_schedule <- function(region, baseline, pars_simulation,
                                     date_end, pars1, booster_eligibility) {
  #browser()
  ## First, the variables that might come from either the simulation
  ## parameters or from the underlying baseline.
  vaccine_uptake <-
    pars_simulation$vaccine_uptake %||% baseline$vaccine_uptake
  if (length(dim(vaccine_uptake)) > 1) {
    ## sircovid::vaccine_schedule_scenario is not yet set up to take
    ## uptake as a matrix, we will just use 1st dose uptake at the moment
    vaccine_uptake <- vaccine_uptake[, 1]
  }

  ## The first clause here should be removed after we have new fits
  ## run after this PR is merged (due to changes in the fitting code).
  ## Note we use `[[` and not `$` to avoid partially matching on
  ## vaccine_eligibility_min_age.  It does look like this won't be an
  ## issue though as we seem to have this available in pars_simulation
  ## atm?
  if (is.null(pars_simulation$vaccine_eligibility_min_age) &&
      !is.null(baseline[["vaccine_eligibility"]])) {
    vaccine_eligibility <- baseline[["vaccine_eligibility"]]
  } else {
    vaccine_eligibility_min_age <-
      pars_simulation$vaccine_eligibility_min_age %||%
      baseline$vaccine_eligibility_min_age
    vaccine_eligibility <-
      sircovid::vaccine_eligibility(vaccine_eligibility_min_age)
  }

  ## This one can't vary in the baseline because it would never have
  ## appeared in the fits.
  doses_future <- pars_simulation$vaccine_daily_doses[[region]]

  ## Assumed never to vary with simulations
  mean_days_between_doses <- baseline$vaccine_mean_days_between_doses

  ## We just need this one to stitch together the full schedule, can't
  ## vary by simulation of course.
  schedule_past <- baseline$vaccine_schedule

  schedule_past$doses <- schedule_past$doses[, -3, ]
  ## These don't vary at all and come from the parameters themselves.
  vaccine_index_dose2 <- pars1$index_dose[[2]]
  vaccine_progression_rate <- pars1$vaccine_progression_rate_base
  N_tot <- pars1$N_tot
  n_strains <- pars1$n_strains

  ## Below here is calculation
  priority_population <- sircovid::vaccine_priority_population(
    region, vaccine_uptake * vaccine_eligibility)

  sircovid::vaccine_schedule_scenario(
    schedule_past = schedule_past,
    doses_future = doses_future,
    end_date = date_end,
    mean_days_between_doses = mean_days_between_doses,
    priority_population = priority_population,
    #boosters_future = pars_simulation$booster_daily_doses[[region]],
    boosters_future = NULL,
    boosters_prepend_zero = FALSE,
    booster_proportion = booster_eligibility)
}


## This all will be worth simplifying, but it looks like it will
## affect the fits too. They will also benefit from moving to using
## the sircovid parameters though as they should always be available
## at that point.

## vaccine_efficacy,
##                                            booster_efficacy, n_strain,
##                                            strain_vaccine_efficacy,
##                                            strain_vaccine_booster_efficacy,
##                                            strain_cross_immunity
simulate_calculate_vaccination <- function(state, index, pars) {
  n_groups <- sircovid:::lancelot_n_groups()
  regions <- dimnames(state)[[3]]
  ## TODO: not clear if we always know (or need to know) where
  ## boosters are in the matrix.
  n_strain <- pars$n_strains # unclear if this should be 2 or 4
  cross_immunity <- pars$cross_immunity # unclear what length we expect

  ## output the cumulative transitions between vaccine strata
  ## by age / vaccine stratum / region / over time
  n_vaccinated <- apply(state[names(index$n_vaccinated), , , , drop = FALSE],
                        c(1, 3, 4), mean)
  n_strata <- nrow(n_vaccinated) / n_groups
  n_vaccinated <-
    mcstate::array_reshape(n_vaccinated, 1L, c(n_groups, n_strata))

  ## output the number recovered in each vaccine stratum / region / over time
  ## R_raw: [age, strain, vaccine, particle, region, time]
  R_raw <- mcstate::array_reshape(
    state[names(index$R), , , , drop = FALSE],
    1L, c(n_groups, n_strain, n_strata))

  ## Take the sum over age
  R <- apply(R_raw, seq(2, 6), sum)
  ## take the mean over the particles
  R <- apply(R, c(1, 2, 4, 5), mean)
  ## R: [strain, vaccine, region, time]

  ## need to allow for imperfect cross-strain immunity
  ## R_strain_1: 100% of strain-level 1, 3, 4 +
  ##   args$cross_protection[2] * strain-level 2
  ## R_strain_2: 100% of strain-level 2, 3, 4 +
  ##   args$cross_protection[1] * strain-level 1
  calc_strain_immunity <- function(strain, R, strain_cross_immunity) {
    apply(R[-strain, , , ], seq(2, 4), sum) +
      R[strain, , , ] * strain_cross_immunity[strain]
  }

  R_strain <- list(strain_1 = calc_strain_immunity(2, R, cross_immunity),
                   strain_2 = calc_strain_immunity(1, R, cross_immunity))


  ## R_strain: [vaccine, region, time]

  # calculate the proportion protected given strain-specific vaccine efficacy
  # and cross-strain immunity
  n_protected_strain_1 <- calculate_n_protected(
    n_vaccinated, R_strain$strain_1, pars, 1)
  dimnames(n_protected_strain_1)[[2]] <- regions

  n_protected_strain_2 <- calculate_n_protected(
    n_vaccinated, R_strain$strain_2, pars, 2)
  dimnames(n_protected_strain_1)[[2]] <- regions

  # Output number of first, second and booster doses

  idx_doses <- c("first_dose" = 1, "second_dose" = 2, "waned_dose" = 3, "booster_dose" = 4)
  doses <- n_vaccinated[, idx_doses, , ]
  dimnames(doses)[2:3] <- list(names(idx_doses), regions)
  doses_inc <- aperm(apply(doses, c(1, 2, 3), diff), c(2, 3, 4, 1))
  doses_inc <- mcstate::array_bind(array(NA, c(dim(doses_inc)[-4], 1)),
                                   doses_inc)
  colnames(doses_inc) <- paste0(colnames(doses), "_inc")

  n_doses <- abind_quiet(doses, doses_inc, along = 2)

  list(n_vaccinated = n_vaccinated,
       n_protected = list(strain_1 = n_protected_strain_1,
                          strain_2 = n_protected_strain_2),
       n_doses = n_doses)
}


## TODO: overlap considerably with calculate_n_protected
## make R strain specific
calculate_n_protected <- function(n_vaccinated, R, pars, strain) {
  vp <- calculate_vaccine_protection(pars, strain)

  # Methodology: calculate incidence of first / second doses,
  # number in each strata in total,
  # number in each strata who are recovered, use these to calculate proportion
  # protected as shown in main fig A / B.
  # to calculate proportion protected we need
  # s = stratum, r = region, t = time
  # let V be number in each vaccination stage V[a, s, r, t]
  # ler R be number recovered R[s, r, t]
  # from top to bottom of figure 1B
  # - number vaccinated sum_sr(V[s, r, t])
  # - number protected after vaccination against severe disease:
  #   > sum_asr(V[a, s, r, t] * eff_severe[a, s])
  # - number protected after vaccination against infection
  #   > sum_asr(V[a, s, r, t] * eff_inf[a, s])
  # - number protected after infection (this is added to all of the above):
  #   > sum_sr(R[s, r, t])
  # - number protected after infection only:
  #   > sum_r(R[1, r, t]


  ## create array of number in each vaccine stratum
  n_strata <- ncol(n_vaccinated)

  # check arrays are conformable
  stopifnot(ncol(vp[[1]]) == n_strata)

  V <- array(0, dim = dim(n_vaccinated))
  V[, n_strata, , ] <- n_vaccinated[, n_strata - 1L, , ]
  for (i in seq(2, n_strata - 1)) {
    V[, i, , ] <- n_vaccinated[, i - 1, , ] - n_vaccinated[, i, , ]
  }

  sum_sr <- function(x) apply(x, c(2, 3), sum)
  sum_asr <- function(x) apply(x, c(3, 4), sum)

  ret <- list(
    ever_vaccinated = sum_sr(n_vaccinated[, 1, , ]),
    protected_against_infection = sum_asr(c(vp$infection) * V),
    protected_against_severe_disease = sum_asr(c(vp$severe_disease) * V),
    protected_against_death = sum_asr(c(vp$death) * V),
    ever_infected = sum_sr(R),
    ever_infected_unvaccinated = R[1, , ]
  )

  aperm(abind_quiet(ret, along = 3), c(3, 1, 2))
}


calculate_vaccine_protection <- function(pars, strain) {
  efficacy_infection <- 1 - pars$rel_susceptibility[, strain, ]
  efficacy_disease <- efficacy_infection + (1 - efficacy_infection) *
    (1 - pars$rel_p_sympt[, strain, ])
  efficacy_severe_disease <- efficacy_disease + (1 - efficacy_disease) *
    (1 - pars$rel_p_hosp_if_sympt[, strain, ])
  efficacy_death <- efficacy_severe_disease + (1 - efficacy_severe_disease) *
    (1 - pars$rel_p_death[, strain, ])

  list(infection = efficacy_infection,
       disease = efficacy_disease,
       severe_disease = efficacy_severe_disease,
       death = efficacy_death)
}


## TODO: this is adapted from spimalot, but here we already have
## transformed parameters so we should just use them.
spim_population <- function(combined, ignore_uk = FALSE, by_age = TRUE,
                            group = 1:19) {
  pop <- vapply(combined$pars[1, ], function(p) p$N_tot[group],
                numeric(length(group)))
  pop_england <- rowSums(pop[, sircovid::regions("england")])
  if (!ignore_uk) {
    pop_uk <- rowSums(pop[, sircovid::regions("all")])
  } else {
    pop_uk <- rep(NA_real_, length(pop_england))
  }

  df <- data.frame(pop, england = pop_england, uk = pop_uk)
  if (!by_age) {
    df <- colSums(df)
  }

  df
}


## rho = beta / rt
rt_to_rho <- function(combined, final_beta_fits, regions) {

  n_particles <- nrow(combined$rt[["Rt_general"]])
  n_regions <- length(regions)
  nms_strain <- dimnames(combined$rt[["Rt_general"]])[[3L]]

  f <- function(rt) {
    array(
      matrix(final_beta_fits, n_particles, n_regions * length(nms_strain)) /
        matrix(
          combined$rt[[rt]][, regions, ], n_particles,
          n_regions * length(nms_strain)
        ),
      c(n_particles, n_regions, 3),
      list(NULL, regions, nms_strain)
    )
  }
  setNames(lapply(c("Rt_general", "eff_Rt_general"), f), c("rt", "rt_eff"))
}

# 2. Function to sample beta
sample_beta <- function(mean, sd, np) {
  sd[sd == 0] <- 1e-14 ## sd must be >0
  sort(distr6::dstr("Lognormal", mean = mean, sd = sd)$rand(np))
}

create_beta_step <- function(key, regions, dt, rt_future, step_end, final_beta_fits,
                            rho, n_beta_add, pars, control_pars, beta) {
  message(sprintf("Computing beta_step for scenario '%s'", key))

  for (region_index in seq_along(regions)) {
    which <- rt_future$scenario == key
    for (i in seq_len(sum(which))) {
      scen_event <- rt_future[which, ][i, ]

      if (scen_event[["trend"]] == "regional") {
        ## CASE 1 - Regional Rt; Scaled; Relative to betas from fits
        # Select region index and propagate uncertainty from fits
        future_beta <- final_beta_fits[, region_index] * scen_event$value
      } else {
        # Mean over regions and particles
        ## FIXME 1 (RS) - Is this how we aggregate regions in the combined? I
        ##                thought it's a weighted Rt.
        ## FIXME 2 (RS)- I don't think this works when including celtic
        ##               nations. We used to only run simulations for England.

        rho_rt <- rho[[scen_event$rt_type]]

        ## CASE 2 - National Rt; Scaled
        if (scen_event[["method"]] == "scaled") {
          # get the mean beta from fits over region and particle
          mean_beta <- mean(final_beta_fits)
          ## CASE 3 - National Rt; Fixed
        } else {
          # get region beta/rt ratio (mean over particles)
          mean_beta <- mean(rho_rt[, region_index, scen_event$strain])
        }

        mean_beta <- mean_beta * scen_event$value

        sigma <- as.numeric(scen_event[["sigma"]]) *
          mean(rho_rt[, , scen_event$strain])

        # We need to add uncertainty as we have taken the mean over particles
        ## TODO (RS) : Do we need to do this? Could we not propagate uncertainty
        ##   over particles like in the regional trend case?
        future_beta <- sample_beta(mean_beta, sigma, nrow(final_beta_fits))
      }

      beta[, region_index, seq(scen_event$step_start, scen_event$step_end)] <- future_beta
    }
  }
  ## All of this moves into spimalot, as a modification of
  ## beta_step, so that the delicate dimension tweaking can be done
  ## in one place.

  ## We need:
  ## * beta_step - the base beta_step computed based on the rt schedule
  ## * rt_seasonality_{amplitide,date_peak}
  ## * rt_schools{modifier,schedule}
  ## * n_beta_add and dt (and/or step_end, step_date, or similar)
  ## * pars or (n_particles, regions)
  ##
  ## which is not that bad.

  ## Back-calculate dates here to deal with any truncating of
  ## step_current:step_end caused by simulating with new data.
  ## Relative distance between us and the peak date (on [0..1])
  ## Simulation step as a fractional date:
  step_date <- seq(to = step_end, length.out = n_beta_add) * dt
  beta_mult_seasonality <- spimalot::spim_beta_mult_seasonality(
    step_date, control_pars$rt_seasonality_date_peak,
    control_pars$rt_seasonality_amplitude)
  ## Replicate appropriately:
  beta_mult_seasonality <-
    array(rep(beta_mult_seasonality, each = prod(dim(beta)[1:2])),
          c(dim(beta)[1:2], n_beta_add))

  ## For schools, we need to reverse-apply this based on school
  ## closure atm in the original fits.  This is likely hard to do
  ## where open/closed is not consistent (but if we just simulate
  ## england that's ok)

  beta_mult_schools <- array(NA_real_, c(dim(pars), length(step_date)))
  for (i in seq_along(regions)) {
    ## Per usual, there's some delicate array expansion here, but
    ## nothing too hairy.  Once this moves into sircovid that'll be
    ## less worrying looking.
    beta_mult_schools_i <- spimalot::spim_beta_mult_schools(
      step_date,
      control_pars$rt_schools_schedule,
      control_pars$rt_schools_modifier,
      regions[[i]]
    )
    beta_mult_schools[, i, ] <- rep(beta_mult_schools_i, each = nrow(pars))
  }

  beta_mult <- beta_mult_seasonality * beta_mult_schools

  ## Then apply into the object
  i <- seq(to = dim(beta)[[3]], length.out = n_beta_add)
  beta[, , i] <- beta[, , i] * beta_mult
  colnames(beta) <- regions

  beta
}

##    a) scale sigma;
##    b) generate Rt sample if needed;
##
## TODO: move into spimalot
compute_beta_step <- function(combined, rt_future, control) {
  ## Start by getting all shared variables

  pars <- combined$pars
  stopifnot(is.list(pars), is.matrix(pars), !is.null(colnames(pars)))

  regions <- colnames(combined$pars)
  scenario <- unique(rt_future$scenario)
  dt <- combined$dt

  step_end <- sircovid::as_sircovid_date(control$date_end) *
              combined$steps_per_day

  beta <- t(vapply(pars, "[[", numeric(length(pars[[1]]$beta_step)),
                  "beta_step"))
  n_beta_add <- step_end - ncol(beta)
  beta <- cbind(beta, matrix(rep(beta[, ncol(beta)], n_beta_add), nrow(beta)))
  beta <- array(beta, c(dim(pars), ncol(beta)))
  final_beta_fits <- beta[, , combined$step]
  colnames(final_beta_fits) <- regions

  ## beta/Rt ratio
  rho <- rt_to_rho(combined, final_beta_fits, regions)

  rt_future <- rt_future %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(step_start = sircovid::sircovid_date(date) / dt,
                  step_end = c(step_start[-1L] - 1L, step_end)) %>%
    dplyr::ungroup()

  setNames(
    lapply(
      scenario,
      create_beta_step,
      regions, dt, rt_future, step_end, final_beta_fits, rho, n_beta_add, pars,
      control$parameters, beta),
    scenario
  )
}

## TODO: maybe move to spimalot
spim_calc_prop_strain_2 <- function(obj) {
  idx1 <- match("n_strain_1_inc", row.names(obj$state))
  idx2 <- match("n_strain_2_inc", row.names(obj$state))
  prop_strain_2 <- obj$state[idx2, , , , drop = FALSE] /
    (obj$state[idx1, , , , drop = FALSE] + obj$state[idx2, , , , drop = FALSE])
  obj$state <- abind_quiet(obj$state, prop_strain_2, along = 1L)
  row.names(obj$state)[nrow(obj$state)] <- "prop_strain_2"
  obj
}

rt_critical_dates <- function(rt_future, combined) {
  rt_dates <- c(rt_future$date, control$parameters$rt_schools_schedule$date)
  sort(sircovid::sircovid_date(unique(rt_dates[rt_dates > combined$date])))
}


rel_severity_strains <- function(vacc_rel_severity_strain1,
                                 vacc_rel_severity_strain2) {
  strain_severity_modifier <- rep(list(list(
    rel_susceptibility = 1,
    rel_p_sympt = 1,
    rel_p_hosp_if_sympt = 1,
    rel_infectivity = 1,
    rel_p_death = 1
  )), 4)
  
  ## modify_severity modifies severity to make it relative to alpha/wildtype according to;
  ## efficacy, efficacy of strain 2 and the strain severity modifier
  ## see `sircovid\R\vaccination.R` for more details
  sircovid::modify_severity(vacc_rel_severity_strain1,
                            vacc_rel_severity_strain2,
                            strain_severity_modifier)
}

tidy_state_one_LOCAL <- function (x, common) 
{
  #browser()
  stopifnot(nrow(common) == 1L)
  res <- list()
  if ("summary_state" %in% names(x)) {
    if (is.null(colnames(x$summary_state))) {
      name_2 <- "particle"
    }
    else {
      name_2 <- "quantile"
    }
    s <- aperm(x$summary_state, c(4, 2, 3, 1))
    dn <- append(dimnames(s), list("all", "all"), 2)
    names(dn) <- c("date", name_2, "group", "vaccine_status", 
                   "region", "state")
    dn$date <- sircovid::sircovid_date_as_date(max(x$date))
    if (is.null(dn[[2]])) {
      dn[[2]] <- seq_len(ncol(s))
    }
    ret <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret$value <- c(s)
    res$summary_state <- ret
  }
  if ("state" %in% names(x)) {
    if (is.null(colnames(x$state))) {
      name_2 <- "particle"
    }
    else {
      name_2 <- "quantile"
    }
    s <- aperm(x$state, c(4, 2, 3, 1))
    dn <- append(dimnames(s), list("all", "all"), 2)
    names(dn) <- c("date", name_2, "group", "vaccine_status", 
                   "region", "state")
    dn$date <- sircovid::sircovid_date_as_date(x$date)
    if (is.null(dn[[2]])) {
      dn[[2]] <- seq_len(ncol(s))
    }
    ret <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret$value <- c(s)
    if (is.null(x$n_protected)) {
      ret_p <- NULL
    }
    else {
      p <- aperm(abind_quiet(x$n_protected, along = 4), 
                 c(3, 2, 4, 1))
      dnp <- set_names(dimnames(p), c("date", "region", 
                                      "strain", "state"))
      dnp$date <- sircovid::sircovid_date_as_date(x$date)
      dnp$quantile <- "mean"
      ret_p <- do.call(expand.grid, dnp)[seq(length(dnp), 
                                             1)]
      ret_p$value <- c(p)
    }
    # if (is.null(x$n_doses)) {
    #   ret_d <- NULL
    # }
    # else {
    #   d <- aperm(x$n_doses, c(4, 1, 3, 2))
    #   dn$group <- c(sircovid:::sircovid_age_bins()$start, 
    #                 "CHW", "CHR")
    #   dn$state <- colnames(x$n_doses)
    #   dn$quantile <- "mean"
    #   ret_d <- do.call(expand.grid, dn)[seq(length(dn), 
    #                                         1)]
    #   ret_d$value <- c(d)
    # }
    ##Removed the above because we're exporting this object in the combined.
    ret_d <- NULL
    
    av <- unlist(lapply(x$state_by_age, aperm, c(4, 1, 2, 
                                                 3)), use.names = FALSE)
    dn$group <- rownames(x$state_by_age[[1]])
    dn$vaccine_status <- colnames(x$state_by_age[[1]])
    dn$state <- names(x$state_by_age)
    ret_av <- do.call(expand.grid, dn)[seq(length(dn), 1)]
    ret_av$value <- av
    res$state <- ret
    res$n_protected <- ret_p
    res$n_doses <- ret_d
    res$state_by_age <- ret_av
    rownames(common) <- NULL
  }
  lapply(res, function(x) cbind(common, x))
}
