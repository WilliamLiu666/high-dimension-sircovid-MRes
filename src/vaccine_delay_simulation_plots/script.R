source("global_util.R")

version_check("sircovid", "0.13.15")
version_check("spimalot", "0.7.11")
##################################################
#FIGURE 1
##################################################
#Made of three sub-plots separately; f1p1, f1p2, f1p3
#f1p1 - cumulative doses over time. Showing three lines; a) all doses
  # b) second doses in the real model fits
  # c) second doses in the 3-week delay counterfactual
#f1p2 - Daily hospital admissions for model fit and counterfactual
#f1p3 - Cumulative hospital admissions for model fit and counterfactual

#Load data
simulation_output <- readRDS('simulation.rds')
fit_trajectories <- readRDS('fit_trajectories.rds')
Real_vaccination <- readRDS('Real_vaccine_schedule.rds')
Counterfactual_vaccination <- readRDS('Counterfactual_vaccine_schedule.rds')
combined <- readRDS('combined.rds')
p_C_and_H <- readRDS('p_C_and_H.rds')
#waning_fit_data <- readRDS("waning_fit_data.rds")

#Prep dataframes for each subplot

f1p1_data <- prep_f1p1_data(Real_vaccination, Counterfactual_vaccination)
f1p2_data <- prep_f1p2_data(simulation_output, fit_trajectories, combined)
f1p3_data <- f1p2_data$Cumulative_data

#PLOTTING

p1 <- plot_f1p1(f1p1_data)
p2 <- plot_f1p2(f1p2_data)
p3 <- plot_f1p3(f1p3_data)

###########


#Plot all three together
f1 <- plot_grid(p1, p2, p3, ncol = 1, align = "v",
                labels = c("A", "B", "C"))


dir.create("figs", FALSE, TRUE)

ggsave(filename = "f1_trajectories.png",
       path = 'figs/', plot = f1,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f1_trajectories.tiff",
       path = 'figs/', plot = f1,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f1_trajectories.pdf",
       path = 'figs/', plot = f1,
       dpi=300, height=11, width=11, units="in")

##################################################
#FIGURE 2
##################################################
#Made of two sub-plots separately; f2p1, f2p2
#f2p1 - vaccine status of entire population over time. 0/1/2/waned dose
#f2p2 - plot showing population average protection due to vaccines over time

#Prepare data
###############
f2p1_data <- prep_f2p1_data(simulation_output, combined)
f2p2_data <- f2p1_data
f2p3_data <- prep_f2p3_data(simulation_output, combined, p_C_and_H)


#Plotting
##########

f2p1 <- plot_f2p1(f2p1_data)
f2p1_collapsed <- plot_f2p1_collapsed(f2p1_data)
f2p2 <- plot_f2p2(f2p2_data)

f2p2_inf_and_hosp <- plot_f2p2_inf_and_hosp(f2p2_data)

f2p3 <- plot_f2p3(f2p3_data)

f2 <- plot_grid(f2p1, f2p2,f2p3,
                   ncol = 1, align = "v", axis = "lr",
                   #rel_heights = c(1,0.6,0.6),
                rel_heights = c(1,1,1),
                   labels = c("A", "B", "C"),
                hjust = 1)

f2_collapsed <- plot_grid(f2p1_collapsed, f2p2, f2p3,
                ncol = 1, align = "v", axis = "lr",
                #rel_heights = c(1,0.6, 0.6),
                rel_heights = c(1,1,1),
                labels = c("A", "B", "C"),
                hjust = 1)


f2_collapsed <- f2_collapsed + theme(plot.margin = unit(c(0,0,0,0.5), "cm")) 

f2 <- f2 + theme(plot.margin = unit(c(0,0,0,0.5), "cm")) 

ggsave(filename = "f2_protection.png",
       path = 'figs/', plot = f2,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f2_protection.tiff",
       path = 'figs/', plot = f2,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f2_protection.pdf",
       path = 'figs/', plot = f2,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f2_collapsed_protection.png",
       path = 'figs/', plot = f2_collapsed,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f2_collapsed_protection.tiff",
       path = 'figs/', plot = f2_collapsed,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f2_collapsed_protection.pdf",
       path = 'figs/', plot = f2_collapsed,
       dpi=300, height=11, width=11, units="in")


##################################################
#FIGURE 3
##################################################
#A facet plot showing daily hospitalisations split by age group


#Prepare data
###############
f3_data <- prep_f3_data(simulation_output, combined)

#Plotting
##########
f3 <- plot_f3(f3_data)
f3_collapsed <- plot_f3_collapsed(f3_data)

f3_SI_stacked <- plot_f3_SI_stacked(f3_data)

f3_age_prop <- plot_f3_SI_age_prop(f3_data)

f3_hospitalisation_distribution <- plot_f3_SI_hosp_dist(f3_data)

ggsave(filename = "f3_hosps.png",
       path = 'figs/', plot = f3,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosps.tiff",
       path = 'figs/', plot = f3,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosps.pdf",
       path = 'figs/', plot = f3,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosps_collapsed.png",
       path = 'figs/', plot = f3_collapsed,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosps_collapsed.tiff",
       path = 'figs/', plot = f3_collapsed,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosps_collapsed.pdf",
       path = 'figs/', plot = f3_collapsed,
       dpi=300, height=7, width=11, units="in")


ggsave(filename = "f3_SI_stacked.png",
       path = 'figs/', plot = f3_SI_stacked,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_SI_proportion.png",
       path = 'figs/', plot = f3_age_prop,
       dpi=300, height=7, width=11, units="in")

ggsave(filename = "f3_hosp_distribution.png",
       path = 'figs/', plot = f3_hospitalisation_distribution,
       dpi=300, height=7, width=11, units="in")



###################
#FIGURE 4 - The waning plot

f4_data_delta <- read.csv("vaccine_efficacy_delta_sensitivities.csv")
f4_data_alpha <- read.csv("vaccine_efficacy_alpha.csv")

f4 <- plot_f4(f4_data_delta, f4_data_alpha)

ggsave(filename = "f4.png",
       path = 'figs/', plot = f4,
       dpi=300, height=7, width=8, units="in")

ggsave(filename = "f4.tiff",
       path = 'figs/', plot = f4,
       dpi=300, height=7, width=8, units="in")

ggsave(filename = "f4.pdf",
       path = 'figs/', plot = f4,
       dpi=300, height=7, width=8, units="in")

###################
#FIGURE 5 - Alpha VEs
f5_alpha <- plot_Alpha(f4_data_alpha)

ggsave(filename = "f5_Alpha_VE.png",
       path = 'figs/', plot = f5_alpha,
       dpi=300, height=7, width=8, units="in")

ggsave(filename = "f5_Alpha_VE.tiff",
       path = 'figs/', plot = f5_alpha,
       dpi=300, height=7, width=8, units="in")

ggsave(filename = "f5_Alpha_VE.pdf",
       path = 'figs/', plot = f5_alpha,
       dpi=300, height=7, width=8, units="in")

###################
#FIGURE 6 - VE SENSITIVITIES
#These plots demonstrate the VEs for our sensitivity scenarios

f6a <- plot_f6a(f4_data_delta)


ggsave(filename = "f6a_sensitivities.png",
       path = 'figs/', plot = f6a,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f6a_sensitivities.tiff",
       path = 'figs/', plot = f6a,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f6a_sensitivities.pdf",
       path = 'figs/', plot = f6a,
       dpi=300, height=11, width=11, units="in")

f6b <- plot_f6b(f4_data_alpha)


ggsave(filename = "f6b_alpha_sensitivities.png",
       path = 'figs/', plot = f6b,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f6b_alpha_sensitivities.tiff",
       path = 'figs/', plot = f6b,
       dpi=300, height=11, width=11, units="in")

ggsave(filename = "f6b_alpha_sensitivities.pdf",
       path = 'figs/', plot = f6b,
       dpi=300, height=11, width=11, units="in")

