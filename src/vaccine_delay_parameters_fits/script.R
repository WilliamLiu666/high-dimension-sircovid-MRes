source("global_util.R")
source("global_vaccine.R")

version_check("sircovid", "0.14.11")
version_check("spimalot", "0.8.20")


date <- "2021-09-13"

## Two epochs after starting with a single strain model (without vaccination)
## * early December 2020: vaccination starts, expand vaccine classes 
##                        but without boosters
## * early March 2021: delta appears, expand strains
epoch_dates <- c("2020-12-07", "2021-03-08")

## The restart date needs to be after the last epoch really, but can
## be quite close.  Once omicron has taken over and the fits have
## settled we will probably keep this 4-8 weeks in the past
restart_date <- "2020-12-08"

## Load all parameters from the last run; creates priors, and updates
## new entries into the proposal matrix as needed.
pars <- load_mcmc_parameters(model_type, assumptions, deterministic,
                             multiregion)

## The baselines are always region-specific
regions <- sircovid::regions("england")
baseline <- lapply(regions, create_baseline,
                   date, model_type, restart_date,
                   epoch_dates, pars$info, assumptions)
names(baseline) <- regions

message("Writing parameters_info.csv")
write_csv(pars$info, "parameters_info.csv")
message("Writing parameters_proposal.csv")
write_csv(pars$proposal, "parameters_proposal.csv")
message("Writing parameters_prior.csv")
write_csv(pars$prior, "parameters_prior.csv")

message("Writing parameters_base.rds")
saveRDS(baseline, "parameters_base.rds")

message("Writing parameters_transform.R")
fs::file_copy("R/transform.R",
              "parameters_transform.R", overwrite = TRUE)
