
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
