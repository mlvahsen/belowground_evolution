# Load libraries
library(tidyverse)
library(multcomp)
library(patchwork)
library(here)
library(rCMEM)

# Read in trait data
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))

## Reviewer 1 - monoculture vs polyculture analysis ####
all_traits %>% 
  mutate(composition = paste(diversity, age, sep = "-")) -> all_traits

# Fit model for root distribution parameter
mod_beta <- lm(beta ~ composition + ic_weight + ln_depth + frame, data = all_traits)

# Check for significant effect of composition
car::Anova(mod_beta)

# See which composition groups are different
pairs(emmeans(mod_beta, "composition", type = "response"))
# Significant differences between mono-ancestral and poly-mix as well as
# poly-ancestral and poly-mix

tibble(beta = summary(emmeans(mod_beta, ~ composition))$emmean,
       lower = summary(emmeans(mod_beta, ~ composition))$lower.CL,
       upper = summary(emmeans(mod_beta, ~composition))$upper.CL,
       composition = summary(emmeans(mod_beta, ~composition))$composition) %>% 
  ggplot(aes(x = reorder(composition, beta), y = beta, shape = composition)) +
  geom_point(size = 4, color = "#e31a1c") +
  geom_jitter(data = all_traits, aes(x = composition, y = beta),
              height = 0, width = 0.1, alpha = 0.3, color = "#e31a1c") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#e31a1c") +
  scale_shape_manual(values = c(16,17,16,8,17)) +
  theme_bw(base_size = 14) +
  ylab(expression(paste("root distribution parameter (", beta, ")"))) +
  xlab("composition (diversity-age)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, vjust =0.7, hjust=0.5)) -> beta_poly_fig

# Repeat the same thing for root-to-shoot ratio
mod_rs <- lm(rs ~ composition + ic_weight + ln_depth + frame, data = all_traits)

# Check for significant effect of composition
car::Anova(mod_rs)

# See which composition groups are different
pairs(emmeans(mod_rs, "composition", type = "response"))
# Significant differences between mono-ancestral and poly-modern as well as
# poly-ancestral and poly-modern

tibble(rs = summary(emmeans(mod_rs, ~ composition))$emmean,
       lower = summary(emmeans(mod_rs, ~ composition))$lower.CL,
       upper = summary(emmeans(mod_rs, ~composition))$upper.CL,
       composition = summary(emmeans(mod_rs, ~composition))$composition) %>% 
  ggplot(aes(x = reorder(composition, rs), y = rs, shape = composition)) +
  geom_point(size = 4, color = "#fb9a99") +
  geom_jitter(data = all_traits, aes(x = composition, y = rs),
              height = 0, width = 0.1, alpha = 0.3, color = "#fb9a99") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#fb9a99") +
  scale_shape_manual(values = c(16,17,16,8,17)) +
  theme_bw(base_size = 14) +
  ylab("root-to-shoot ratio") +
  xlab("composition (diversity-age)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, vjust =0.7, hjust=0.5))-> rs_poly_fig

# Bring graphs together and save
png("~/Desktop/poly_fig_revisions.png", height = 4, width = 9, res = 300, units = "in")
rs_poly_fig + beta_poly_fig
dev.off()

## Reviewer 3 - issues with not varying Z_min and Z_max by genotype ####

# Address this in 3 ways: (1) carbon sequestration and marsh accretion are not
# particularly sensitive to Z_min and Z_max, (2) this is particularly true at
# conservative estimates of sea-level rise, and (3) genotypes don't vary that
# much in Z_max and Z_min from other experiments

# (1) The model is not very sensitive to Z_min and Z_max

# We first need to update our tidal data which was incorrect in the original
# submission version of the paper

# Get mean annual tidal data
tides <- read_csv(here("supp_data", "tides_2018.csv"))
mean(tides$`MSL (m)`) -> msl
mean(tides$`MHW (m)`) -> mhw


# Blue genes estimates for Z_min and Z_max (in m)
zMax_for_sim <- 0.5867284
zMin_for_sim <- -0.001624103

agb_cohort_forMEM <- c(0.07794736, 0.08732748, 0.07520756, 0.08458767)
root_shoot_cohort_forMEM <- c(1.789741, 1.541963, 1.649520, 1.401742)
rooting_depth <- c(29.62963, 23.67367, 25.47548, 20.87087)

# Set up storage
run_store_cohort <- matrix(NA, nrow = 4, ncol = 100)
carbon_store_cohort <- matrix(NA, nrow = 4, ncol = 100)

for (i in 1:4){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = agb_cohort_forMEM[i], 
                                 zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = root_shoot_cohort_forMEM[i],
                                 rootTurnover=0.55, rootDepthMax=rooting_depth[i], omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_cohort[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_cohort[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rates_cohort <- (run_store_cohort[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate_cohort <- (carbon_store_cohort[,80] - carbon_store_cohort[,1]) / 80

# Create a data frame to hold all of that information
tibble(location = c("corn", "sellman", "corn", "sellman"),
       age = c("ancestral", "ancestral", "modern", "modern"),
       acc_v = avg_accretion_rates_cohort * 10,
       # unit conversion for carbon accumulation 
       acc_C = avg_C_accum_rate_cohort * 1e-6 / 1e-8) -> cohort_summary

# Now run sensitivity analysis with AGB parameters: 10% increase and decrease in
# each parameter

tibble(zMax_sens = c(zMax_for_sim*0.9, zMax_for_sim, zMax_for_sim*1.1),
       zMin_sens = c(zMin_for_sim*0.9, zMin_for_sim, zMin_for_sim*1.1),
       bMax_sens = c(mean(agb_cohort_forMEM)*0.9, mean(agb_cohort_forMEM),
                     mean(agb_cohort_forMEM)*1.1)) -> sensitivity_params

# Run 1 - Change bMax and hold all parameters constant

# Set up storage
run_store_bMax <- matrix(NA, nrow = 3, ncol = 100)
carbon_store_bMax <- matrix(NA, nrow = 3, ncol = 100)

for (i in 1:3){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = sensitivity_params$bMax_sens[i], 
                                 zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                                 rootTurnover=0.55, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_bMax[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_bMax[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rates_bMax <- (run_store_bMax[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate_bMax <- (carbon_store_bMax[,80] - carbon_store_bMax[,1]) / 80

# Create a data frame to hold all of that information
tibble(percentile = c("90%", "100%", "110%"),
       `peak aboveground biomass (g/cm2)` = sensitivity_params$bMax_sens,
       `vert. accretion rate` = avg_accretion_rates_bMax * 10,
       # unit conversion for carbon accumulation 
       `C accum. rate` = avg_C_accum_rate_bMax * 1e-6 / 1e-8) -> bMax_summary

# Run 2 - Change Z_min and hold all others constant

# Set up storage
run_store_zMin <- matrix(NA, nrow = 3, ncol = 100)
carbon_store_zMin <- matrix(NA, nrow = 3, ncol = 100)

for (i in 1:3){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = sensitivity_params$bMax_sens[2], 
                                 zVegMin=sensitivity_params$zMin_sens[i]*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                                 rootTurnover=0.55, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_zMin[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_zMin[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rates_zMin <- (run_store_zMin[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate_zMin <- (carbon_store_zMin[,80] - carbon_store_zMin[,1]) / 80

# Create a data frame to hold all of that information
tibble(percentile = c("90%", "100%", "110%"),
       zMin = sensitivity_params$zMin_sens*100,
       acc_v = avg_accretion_rates_zMin * 10,
       # unit conversion for carbon accumulation 
       acc_C = avg_C_accum_rate_zMin * 1e-6 / 1e-8) -> zMin_summary

# Run 3 - Change Z_max and hold all others constant

# Set up storage
run_store_zMax <- matrix(NA, nrow = 3, ncol = 100)
carbon_store_zMax <- matrix(NA, nrow = 3, ncol = 100)

for (i in 1:3){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = sensitivity_params$bMax_sens[2], 
                                 zVegMin=zMin_for_sim*100, zVegMax=sensitivity_params$zMax_sens[i]*100, zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                                 rootTurnover=0.55, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_zMax[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_zMax[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rates_zMax <- (run_store_zMax[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate_zMax <- (carbon_store_zMax[,80] - carbon_store_zMax[,1]) / 80

# Create a data frame to hold all of that information
tibble(percentile = c("90%", "100%", "110%"),
       zMax = sensitivity_params$zMax_sens*100,
       acc_v = avg_accretion_rates_zMax * 10,
       # unit conversion for carbon accumulation 
       acc_C = avg_C_accum_rate_zMax * 1e-6 / 1e-8) -> zMax_summary

# Run 4 - 10% increase/decrease in elevation range?
(zMax_for_sim*100 - zMin_for_sim*100) * 1.1 -> inc_range
# Calculate the amount to add or subtract from zMin or zMax
(inc_range - (zMax_for_sim*100 - zMin_for_sim*100))/2 -> inc_range_add

# Same for decrease
(zMax_for_sim*100 - zMin_for_sim*100) * 0.9 -> dec_range
# Calculate the amount to add or subtract from zMin or zMax
(dec_range - (zMax_for_sim*100 - zMin_for_sim*100))/2 -> dec_range_add
# This is the same value as it is for increase!

tibble(zMin = c(zMin_for_sim*100 + inc_range_add, zMin_for_sim*100, zMin_for_sim*100 - inc_range_add),
       zMax = c(zMax_for_sim*100 - inc_range_add, zMax_for_sim*100, zMax_for_sim*100 + inc_range_add)) -> sens_range

# Set up storage
run_store_range <- matrix(NA, nrow = 3, ncol = 100)
carbon_store_range <- matrix(NA, nrow = 3, ncol = 100)

for (i in 1:3){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = sensitivity_params$bMax_sens[2], 
                                 zVegMin = sens_range$zMin[i], zVegMax=sens_range$zMax[i], zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                                 rootTurnover=0.55, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_range[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_range[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rates_range <- (run_store_range[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate_range <- (carbon_store_range[,80] - carbon_store_range[,1]) / 80

# Create a data frame to hold all of that information
tibble(percentile = c("90%", "100%", "110%"),
       `range (cm)` = c(dec_range, zMax_for_sim*100 - zMin_for_sim*100, inc_range),
       `vert. accretion rate` = avg_accretion_rates_range * 10,
       # unit conversion for carbon accumulation 
       `C accum. rate` = avg_C_accum_rate_range * 1e-6 / 1e-8) -> range_summary

# Plot this up for presenting
quadraticABC <- function(zVegMin, bMax, zVegMax){
  zVegPeak <- (zVegMax - zVegMin)/2
  a <- -((-zVegMin * bMax - zVegMax * bMax) / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  b <- -(bMax / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  c <- (zVegMin * zVegMax * bMax) / ((zVegMin - zVegPeak) * (zVegMax - zVegPeak))
  return(c(a,b,c))
}

quadraticABC(sens_range$zMin[1], as.numeric(bMax_for_sim), sens_range$zMax[1]) -> quad_params1
quadraticABC(sens_range$zMin[2], as.numeric(bMax_for_sim), sens_range$zMax[2]) -> quad_params2
quadraticABC(sens_range$zMin[3], as.numeric(bMax_for_sim), sens_range$zMax[3]) -> quad_params3

tibble(z = seq(-5,65,0.5),
       inc_range10 = quad_params1[1]*z + quad_params1[2]*z^2 + quad_params1[3],
       normal_range = quad_params2[1]*z + quad_params2[2]*z^2 + quad_params2[3],
       dec_range10 = quad_params3[1]*z + quad_params3[2]*z^2 + quad_params3[3]) %>% 
  gather(key = simulation, value = agb, inc_range10:dec_range10) %>% 
  mutate(agb = ifelse(agb < 0, 0, agb)) %>% 
  ggplot(aes(x = z, y = agb, color = simulation)) +
  geom_line(size = 1.2) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  ylim(0, 0.09) -> range_plot

# Repeat the same for bMax sensitivity
quadraticABC(zMin_for_sim*100, sensitivity_params$bMax_sens[1], zMax_for_sim*100) -> quad_params4
quadraticABC(zMin_for_sim*100, sensitivity_params$bMax_sens[2], zMax_for_sim*100) -> quad_params5
quadraticABC(zMin_for_sim*100, sensitivity_params$bMax_sens[3], zMax_for_sim*100) -> quad_params6

tibble(z = seq(-5,65,0.5),
       bMax90 = quad_params4[1]*z + quad_params4[2]*z^2 + quad_params4[3],
       normal = quad_params5[1]*z + quad_params5[2]*z^2 + quad_params5[3],
       bMax110 = quad_params6[1]*z + quad_params6[2]*z^2 + quad_params6[3]) %>% 
  gather(key = simulation, value = agb, bMax90:bMax110) %>% 
  mutate(agb = ifelse(agb < 0, 0, agb)) %>% 
  ggplot(aes(x = z, y = agb, color = simulation)) +
  geom_line(size = 1.2) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Dark2") +
  ylim(0, 0.09) -> bMax_plot

library(patchwork)
png("~/Desktop/sensitivity.png", height = 3.5, width = 8.5, res = 300, units = "in")
bMax_plot + range_plot
dev.off()

bMax_summary
range_summary
