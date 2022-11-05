## TODO: Someone clean this up, is is needlessly convoluted for what
## looks like a fairly straightforward bit of data manipulation.
## Consider the starting point and ending point and see if it can be
## reduced to a few understandable lines.
##
## TODO: at the moment this assumes there is just a single scenario,
## so will need lots of work before we add more.
FIXME_VACCINE_SCHEDULE <- function(future_doses, regions) {
  message("Fed up with seeing this message? Please fix this function!")

  regional_pop <- vapply(regions, function(r)
    sum(sircovid:::sircovid_population(r)),
    numeric(1))
  regional_prop_pop <- regional_pop / sum(regional_pop)

  future_doses <- future_doses %>%
    tidyr::pivot_longer(!week_start, names_to = "type_doses",
                        values_to = "doses_per_week") %>%
    split(., .$type_doses)
  future_doses_by_region <-
    lapply(future_doses, split_doses_by_region, regional_prop_pop)

  daily_doses <- lapply(regions, function (region)
    list(central = get_daily_doses(
           region, date = date,
           future_doses_by_region = future_doses_by_region$doses_per_week)))
  names(daily_doses) <- regions
  daily_doses <- switch_levels(daily_doses)

  booster_daily_doses <- lapply(regions, function (region)
    list(central = get_daily_doses(
           region, date = date,
           future_doses_by_region = future_doses_by_region$boosters_per_week)))
  names(booster_daily_doses) <- regions
  booster_daily_doses <- switch_levels(booster_daily_doses)

  list(vaccine_daily_doses = daily_doses,
       booster_daily_doses = booster_daily_doses)
}


get_daily_doses <- function(region, date, future_doses_by_region = NULL) {
  ## TODO: this all seems like unused heritage code
  ## check with Raphael to delete
  date_start <- as.character(date)
  date_start_plus_two <- as.character(date + 2)
  date_start_plus_three <- as.character(date + 3)

  future <- round(future_doses_by_region[[region]] / 7)
  future_date <- future_doses_by_region$week_start
  if (any(diff(future_date) <= 0)) {
    stop("future_date$week_start must be strictly increasing")
  }
  min_future_date <- date + 2
  ret <- check_future_dates(future, future_date, min_future_date)

  ## All dates must be increasing; just checks the logic above
  if (length(ret) > 1) {
    stopifnot(all(diff(as.Date(names(ret))) > 0))
  }

  ret
}


check_future_dates <- function(future, future_date, min_future_date) {
  if (any(future_date < min_future_date)) {
    i <- max(which(future_date <= min_future_date))
    future <- future[i:length(future)]
    future_date <- future_date[i:length(future_date)]
    future_date[[1]] <- min_future_date
  }
  names(future) <- as.character(future_date)
  future
}


split_doses_by_region <- function(future_doses,
                                  regional_prop_pop) {
  week_start <- as.Date(future_doses$week_start, "%d-%B-%y")
  future_doses$week_start <- week_start

  message("Using following weekly vaccine roll-out schedule(s): ")
  msg <- apply(future_doses, 1, paste, collapse= ": ")
  message(paste(msg, collapse = "\n"))

  doses_per_week_cols <- grep("doses_per_week", names(future_doses), value = TRUE)

  if (length(doses_per_week_cols) > 1) {
    doses_per_week <- lapply(doses_per_week_cols, function(e) future_doses[[e]])
    names(doses_per_week) <- gsub("doses_per_week_", "", doses_per_week_cols)
  } else {
    doses_per_week <- list(main = future_doses$doses_per_week)
  }

  ### split future_doses by region
  future_doses_by_region <- lapply(doses_per_week, function(e) {
    x <- as.data.frame(matrix(
      vapply(regional_prop_pop, function(r) e * r, numeric(length(e))),
      nrow = length(e)))
    names(x) <- names(regional_prop_pop)
    x$week_start <- future_doses$week_start
    x
  })

  if (length(future_doses_by_region) > 1) {
    future_doses_by_region
  } else {
    future_doses_by_region$main
  }
}


check_schools_schedule <- function(schools_schedule, region, target_date) {

  test_nation <-
    ifelse(region %in% sircovid::regions("england"), "england", region)

  test_schedule <- schools_schedule %>%
    dplyr::filter(nation == test_nation)


  test_date <- as.Date(with(test_schedule,
                            paste(year, month, day, sep="-")), "%Y-%m-%d")

  test_date <- tail(which(target_date >= test_date), 1)

  string_out <- paste0("schools_", test_schedule[test_date, "schools"])
}


vacc_efficacy_delta_omicron <- function(data_vaccination,
                                        ve_delta, ve_omicron) {
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
  # 100% PF in 0-40
  prop_pfizer[sircovid:::sircovid_age_bins()$end < 40] <- 1

  ## Get efficacy for Delta
  efficacy_each_vaccine_class <- calculate_average_vacc_efficacy(ve_delta, prop_pfizer)
  vacc_efficacy_delta <- lapply(efficacy_each_vaccine_class, function(e)
    get_vaccine_conditional_prob(e$death,
                                 e$severe_disease,
                                 e$disease,
                                 e$infection,
                                 e$transmission))$central

  efficacy_each_vaccine_class_omicron <-
    calculate_average_vacc_efficacy(ve_omicron, prop_pfizer)
  vacc_efficacy_omicron <- lapply(
    efficacy_each_vaccine_class_omicron, function(e)
      get_vaccine_conditional_prob(e$death,
                                   e$severe_disease, e$disease,
                                   e$infection, e$transmission))

    strain_severity_modifier_omicron <- rep(list(list(
      rel_susceptibility = 1,
      rel_p_sympt = 1,
      rel_p_hosp_if_sympt = 1,
      rel_infectivity = 1,
      rel_p_death = 1
    )), 4)

    ## Apply a 1.85 multiplier to both Delta and Omicron (and metastrains)
    ## relative to Alpha
    strain_severity_modifier_omicron[[1]]$rel_p_hosp_if_sympt <- 1.85
    strain_severity_modifier_omicron[[2]]$rel_p_hosp_if_sympt <- 1.85
    strain_severity_modifier_omicron[[3]]$rel_p_hosp_if_sympt <- 1.85
    strain_severity_modifier_omicron[[4]]$rel_p_hosp_if_sympt <- 1.85

    strain_severity_modifier_omicron <-
      rep(list(strain_severity_modifier_omicron), length(vacc_efficacy_omicron))
    names(strain_severity_modifier_omicron) <- names(vacc_efficacy_omicron)

    for (nm in names(vacc_efficacy_omicron)) {
      strain_severity_modifier_omicron[[nm]][[3]]$rel_p_hosp_if_sympt <-
        strain_severity_modifier_omicron[[nm]][[3]]$rel_p_hosp_if_sympt *
        (1 - ve_omicron[[nm]][ve_omicron$vaccine == "PF" &
                                ve_omicron$type == "severe_disease" &
                                ve_omicron$dose == 2]) /
        (1 - ve_omicron[[nm]][ve_omicron$vaccine == "PF" &
                                ve_omicron$type == "disease" &
                                ve_omicron$dose == 2])
    }

    rel_efficacy <- lapply(names(vacc_efficacy_omicron), function(x)
      sircovid::modify_severity(vacc_efficacy_delta, vacc_efficacy_omicron[[x]],
                                strain_severity_modifier_omicron[[x]]))
    names(rel_efficacy) <- names(vacc_efficacy_omicron)

    rel_efficacy
}

## checked required assertions and set starting dates and standard deviations
##  (if not given)
read_rt_scenarios <- function(path, date, sigma) {

  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d[d$date == "today", "date"] <- as.character(date)
  d$date <- as.Date(d$date)
  if (any(d$date < date)) {
    stop(sprintf("'date' must be >=%s", date))
  }
  ## set sigma from other parameters if not given
  d[is.na(d$sigma), "sigma"] <- sigma

  ## strictly speaking this assertion isn't required however adding it should
  ##  catch the problem of accidentally omitting scenarios where R doesn't
  ##  change (i.e. MTPs)
  start_len <- d %>%
    dplyr::group_by(scenario) %>%
    dplyr::filter(date == min(date), event != "start") %>%
    nrow()
  if (start_len > 0) {
    stop("First event of all scenarios must be called 'start'")
  }

  ## TODO - ADD CHECK THAT ALL DATES FOR ONE SCENARIO ARE DIFFERENT
  n_duplicated_dates <- d %>%
    dplyr::group_by(scenario) %>%
    dplyr::count(date) %>%
    dplyr::filter(n > 1) %>%
    nrow()
  if (n_duplicated_dates > 0) {
    stop("Dates within scenarios must be unique")
  }

  if (!all(d$method %in% c("scaled", "fixed"))) {
    stop("'method' must be 'scaled' or 'fixed'")
  }

  if (!all(d$trend %in% c("regional", "national"))) {
    stop("'trend' must be 'regional' or 'national'")
  }

  if (!all(d$strain %in% c("both", "strain_1", "strain_2"))) {
    stop("'strain' must be one of 'both', 'strain_1', 'strain_2'")
  }

  if (!all(d$rt_type %in% c("rt", "rt_eff"))) {
    stop("'rt_type' must be 'rt' or 'rt_eff'")
  }

  if (any(d$value <= 0)) {
    stop("'value' must be >0")
  }

  if (any(d$sigma < 0)) {
    stop("'sigma' must be >=0")
  }

  d
}
