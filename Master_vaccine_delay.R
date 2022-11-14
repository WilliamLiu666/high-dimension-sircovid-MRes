
## ---------------------------
root_dir <- paste0(orderly::orderly_config()$root, "/src/")
## ---------------------------

#------------------------------------------------------------------------------
#DATA SETUP TASK
fits_data <- orderly::orderly_run(
  "vaccine_delay_fits_data"
)

orderly::orderly_commit(fits_data)
#------------------------------------------------------------

## Parameters fits task sets up before fitting for each region
pars_bb <- orderly::orderly_run(
  "vaccine_delay_parameters_fits",
  parameters = list(deterministic = FALSE),
  use_draft = "newer")
orderly::orderly_commit(pars_bb)


### ---------------------------------------------------------------------------
##--------------------
## SHORT RUNS
# Set as TRUE for short runs to investigate and de-bug
# Set as FALSE for far longer runs and better fits
short_run <- FALSE

##DETERMINISTIC
# Model can be run deterministically for quicker de-bugging process or fit tuning
# Note that simulation tasks downstream require stochastic (deterministic <- FALSE)
deterministic <- FALSE
##--------------------
regions <- sircovid::regions("england")

#-------------
#RUN
bb_fits_ids <- lapply(X = regions,
                 FUN = function(x) {
                   orderly::orderly_run('vaccine_delay_fits',
                                        parameters = list(region = x,
                                                          short_run = short_run,
                                                          deterministic = deterministic),
                                        use_draft = "newer")})

lapply(X= bb_fits_ids,
       FUN = function(x){
         orderly::orderly_commit(x)
       })


# combine regions
combined_bb <- orderly::orderly_run('vaccine_delay_fits_combined',
                                    parameters = list(short_run = short_run,
                                                      deterministic = deterministic
                                                      ),
                                    use_draft = "newer")

orderly::orderly_commit(combined_bb)
