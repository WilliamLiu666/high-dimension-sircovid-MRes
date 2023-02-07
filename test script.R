#This script is purely used for testing on the student's pc.

setwd("C:/Users/WILL LIU/Desktop/Data Science/high-dimension-sircovid-MRes")

orderly::orderly_run(
  "vaccine_delay_fits_data"
)

orderly::orderly_run(
  "vaccine_delay_parameters_fits",
  parameters = list(deterministic = TRUE),
  use_draft = "newer")

orderly::orderly_develop_start('vaccine_delay_fits_gradient',
                     parameters = list(region = "london",
                                       short_run = TRUE,
                                       deterministic = TRUE),
                     use_draft = "newer")

orderly::orderly_run('vaccine_delay_fits_gradient',
                               parameters = list(region = "london",
                                                 short_run = TRUE,
                                                 deterministic = TRUE),
                               use_draft = "newer")

orderly::orderly_develop_start('vaccine_delay_fits_gradient_comp',
                               parameters = list(region = "london",
                                                 short_run = TRUE,
                                                 deterministic = TRUE),
                               use_draft = "newer")
