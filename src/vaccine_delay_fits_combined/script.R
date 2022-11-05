source("global_util.R")

version_check("sircovid", "0.13.15")
version_check("spimalot", "0.7.11")

date <- "2021-09-13"
#date <- "2020-12-08"

sircovid_model <- "lancelot"
model_type <- "BB"

dat <- spimalot::spim_combined_load("regional_results", "england")



dir.create("outputs", FALSE, TRUE)
dir.create("figs", FALSE, TRUE)
dir.create("figs_by_age", FALSE, TRUE)
dir.create("spim_view", FALSE, TRUE)

#Paste in restart state and inflate to match end-of-fits epoch
onward_2020_12_08 <- combined_onward_restart(dat$onward, "2020-12-08",
                                             "regional_results",
                                             "england")


saveRDS(dat$data, "outputs/aggregated_data.rds")

saveRDS(dat$rt$england, "regional_results/Rt_england.rds")

#Instead of saving the original onward state from end of fits, we instead output
#the state from our restart date.
#saveRDS(dat$onward, "outputs/combined.rds")
saveRDS(onward_2020_12_08, "outputs/combined.rds")

spimalot::spim_pars_pmcmc_save(dat$parameters, "outputs/parameters")

write_png("figs/forest_plot.png",
          width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "non_betas"))

write_png("figs/forest_plot_betas.png", width = 2400, height = 1600, res = 200,
          spim_plot_forest(dat, plot_type = "betas"))

write_png("figs/data_fits_regional.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths_hosp", "deaths_carehomes", "deaths_comm", "icu",
              "general", "hosp", "all_admission"), age_band = "all",
            with_forecast = FALSE, add_betas = FALSE))


write_png("figs/serology_euroimmun.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(dat, sircovid::regions("england"), 1, 40))

write_png("figs/serology_roche_n.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_serology(dat, sircovid::regions("england"), 2, 40))


pillar2_age_bands <-
  c("over25", "under15", "15_24", "25_49", "50_64", "65_79", "80_plus")
if (model_type == "BB") {
  write_png("figs/pillar2_all_ages.png", width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_pillar2_positivity(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date("2020-05-15"), ymax = 50))
  for (i in pillar2_age_bands) {
    if (i == "over25") {
      fig_name <- "figs/pillar2_over25.png"
    } else if (i == "under15") {
      fig_name <- "figs_by_age/pillar2_0_14.png"
    } else {
      fig_name <- paste0("figs_by_age/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_positivity(
                dat, sircovid::regions("england"), i,
                date_min = as.Date("2020-05-15"), ymax = 50))
  }
} else if (model_type == "NB") {
  write_png("figs/pillar2_all_ages.png", width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_pillar2_cases(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date("2020-05-15")))
  for (i in pillar2_age_bands) {
    if (i == "over25") {
      fig_name <- "figs/pillar2_over25.png"
    } else if (i == "under15") {
      fig_name <- "figs_by_age/pillar2_0_14.png"
    } else {
      fig_name <- paste0("figs_by_age/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_cases(
                dat, sircovid::regions("england"), i,
                date_min = as.Date("2020-05-15")))
  }
}

write_png("figs/react.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"), date_min = as.Date("2020-05-15"),
            ymax = 10))

write_png("figs/variant_Alpha_Delta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_variant(
            dat, sircovid::regions("england"), "Delta",
            date_min = as.Date("2021-03-01"),
            date_max = as.Date("2021-08-15")))

write_png("figs/incidence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england")))

write_png("figs/incidence_per_1000.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_incidence(
            dat, c(sircovid::regions("england"), "england"), per_1000 = TRUE))

write_png("figs/Rt_all.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "Rt_all"))

write_png("figs/Rt_eff_all.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "eff_Rt_all"))

write_png("figs/Rt_eff_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"),
            "eff_Rt_general"))

write_png("figs/Rt_general.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "Rt_general"))

write_png("figs/beta.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_Rt(
            dat, c(sircovid::regions("england"), "england"), "beta"))

## plotting admissions demography
write_png("figs/admissions_demo.png", width = 2400, height = 1000, res = 200,
          spimalot::spim_plot_admissions_by_age(dat, "england"))


## add (zoomed in) plots of SPI-M-relevant trajectories
write_png("spim_view/regions.png", width = 2400 / 5 * 7, height = 1800,
          res = 200,
          spimalot::spim_plot_trajectories(
            dat, sircovid::regions("england"),
            c("deaths", "deaths_hosp", "icu", "general",
              "hosp", "all_admission"),
            date_min = as.Date(dat$info$date) - 75, age_band = "all",
            with_forecast = FALSE, add_betas = TRUE))

if (model_type == "BB") {
  write_png("spim_view/pillar2_all_ages.png", width = 2400, height = 1200,
            res = 200,
            spimalot::spim_plot_pillar2_positivity(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date(dat$info$date) - 75,
              ymax = 50, add_betas = TRUE))
  for (i in pillar2_age_bands) {
    if (i == "under15") {
      fig_name <- "spim_view/pillar2_0_14.png"
    } else {
      fig_name <- paste0("spim_view/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_positivity(
                dat, sircovid::regions("england"), i,
                date_min = as.Date(dat$info$date) - 75,
                ymax = 50, add_betas = TRUE))
  }
} else if (model_type == "NB") {
  write_png("spim_view/pillar2_all_ages.png", width = 2400, height = 1200,
            res = 200,
            spimalot::spim_plot_pillar2_cases(
              dat, sircovid::regions("england"), "all",
              date_min = as.Date(dat$info$date) - 75,
              add_betas = TRUE))
  for (i in pillar2_age_bands) {
    if (i == "under15") {
      fig_name <- "spim_view/pillar2_0_14.png"
    } else {
      fig_name <- paste0("spim_view/pillar2_", i, ".png")  
    }
    write_png(fig_name, width = 2400, height = 1200, res = 200,
              spimalot::spim_plot_pillar2_cases(
                dat, sircovid::regions("england"), i,
                date_min = as.Date(dat$info$date) - 75,
                add_betas = TRUE))
  }
}

write_png("spim_view/prevalence.png", width = 2400, height = 1200, res = 200,
          spimalot::spim_plot_react(
            dat, sircovid::regions("england"),
            date_min = as.Date(dat$info$date) - 75,
            ymax = 10, add_betas = TRUE))

## Plot outputs by age

deaths_hosp_age_bands <- c("0_49", "50_54", "55_59", "60_64", "65_69", "70_74",
                           "75_79", "80_plus")
for (i in deaths_hosp_age_bands) {
  fig_name <- paste0("figs_by_age/deaths_hosp_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, sircovid::regions("england"), "deaths_hosp", age_band = i,
              with_forecast = FALSE, add_betas = FALSE))
  
  fig_name <- paste0("spim_view/deaths_hosp_", i, ".png")
  write_png(fig_name, width = 2400, height = 1200, res = 200,
            spimalot::spim_plot_trajectories_by_age(
              dat, sircovid::regions("england"), "deaths_hosp", age_band = i,
              date_min = as.Date(dat$info$date) - 75,
              with_forecast = FALSE, add_betas = FALSE))
}
