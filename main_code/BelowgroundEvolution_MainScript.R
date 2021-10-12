# Vahsen et al. Belowground Evolution main script

# This script runs the main analyses and sources additional functions within the
# main_code folder

## Load libraries ####
library(tidyverse); library(ggmcmc); library(rjags); library(lme4); library(here)


## Read in data ####

# Source in initial conditions (genotypes and wet weights)
source(here("supp_code", "GCREW_InitialConditions.R"))
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))

# Set all factors
all_traits %>% 
  mutate(pot = factor(pot),
         frame = factor(frame)) -> all_traits

# Subset out monocultures for monoculture trait analysis
mono_traits <- all_traits %>% 
  filter(diversity == "mono")


## Source JAGS functions to fit (generalized) linear mixed models ####
source("main_code/jags_fxns.R")


## Create a model template that JAGS will use the model matrix for ####
model_template_mono <- lmer(agb ~ ln_depth + ic_weight + frame +
                         (1|genotype), data = mono_traits)


## Set MCMC specifications ####
n.iter <- 10000
n.adapt <- 1000
thin <- 3


## Fit all basic trait models ####

## aboveground biomass
agb_out <- run_models(mono_traits, "agb", model_template_mono, diag_plot = F)
## stem density
density_out <- run_models(mono_traits, "density", model_template_mono, diag_plot = F)
# Fails to set monitor for sigma.res because there is no residual variation in
# Poisson regression (aka this warning is ok and expected)
## stem height
height_out <- run_models(mono_traits, "mean_tot_height", model_template_mono, diag_plot = F)
## stem width
width_out <- run_models(mono_traits, "mean_mid_width", model_template_mono, diag_plot = F)
## belowground biomass
bgb_out <- run_models(mono_traits, "bgb", model_template_mono, diag_plot = F)
## root-to-shoot ratio
rs_out <- run_models(mono_traits, "rs", model_template_mono, diag_plot = F)
## beta (belowground biomass distribution parameter)
beta_out <- run_models(mono_traits, "beta", model_template_mono, diag_plot = F)

# Save as R data file objects to use later for plotting and subsequent analyses
saveRDS(agb_out, here("outputs/monoculture_models", "agb_monomodel.rds"))
saveRDS(rs_out, here("outputs/monoculture_models", "rs_monomodel.rds"))
saveRDS(bgb_out, here("outputs/monoculture_models", "bgb_monomodel.rds"))
saveRDS(width_out, here("outputs/monoculture_models", "width_monomodel.rds"))
saveRDS(height_out, here("outputs/monoculture_models", "height_monomodel.rds"))
saveRDS(density_out, here("outputs/monoculture_models", "density_monomodel.rds"))
saveRDS(beta_out, here("outputs/monoculture_models", "beta_monomodel.rds"))


## Fit all trait models for age + provenance additive effects ####

# Subset data to just Corn Island and Sellman Creek because those provenances
# have multiple reps within each age cohort (only drops 6 reps)
mono_traits %>% 
  filter(location %in% c("corn", "sellman")) -> cs_traits

# Previous models fits showed non-significant age by provenance interactions for
# all trait models. So here we fit all models with additive effects of age and
# provenance 

# Create new model template
model_template_cs <- lmer(agb ~ ln_depth + ic_weight + frame +
                            location + age + (1|genotype), data = cs_traits)

## aboveground biomass
agb_cs_out <- run_models(cs_traits, "agb", model_template_cs, diag_plot = F)
## stem density
density_cs_out <- run_models(cs_traits, "density", model_template_cs, diag_plot = F)
# Fails to set monitor for sigma.res because there is no residual variation in
# Poisson regression (aka this warning is ok and expected)
## stem height
height_cs_out <- run_models(cs_traits, "mean_tot_height", model_template_cs, diag_plot = F)
## stem width
width_cs_out <- run_models(cs_traits, "mean_mid_width", model_template_cs, diag_plot = F)
## belowground biomass
bgb_cs_out <- run_models(cs_traits, "bgb", model_template_cs, diag_plot = F)
## root-to-shoot ratio
rs_cs_out <- run_models(cs_traits, "rs", model_template_cs, diag_plot = F)
## beta (belowground biomass distribution parameter)
beta_cs_out <- run_models(cs_traits, "beta", model_template_cs, diag_plot = F)

# Save as R data file objects to use later for plotting and subsequent analyses
saveRDS(agb_cs_out, here("outputs/corn_sellman_models/", "agb_csmodel.rds"))
saveRDS(rs_cs_out, here("outputs/corn_sellman_models/", "rs_csmodel.rds"))
saveRDS(beta_cs_out, here("outputs/corn_sellman_models/", "beta_csmodel.rds"))
saveRDS(width_cs_out, here("outputs/corn_sellman_models/", "width_csmodel.rds"))
saveRDS(height_cs_out, here("outputs/corn_sellman_models/", "height_csmodel.rds"))
saveRDS(bgb_cs_out, here("outputs/corn_sellman_models/", "bgb_csmodel.rds"))
saveRDS(density_cs_out, here("outputs/corn_sellman_models/", "density_csmodel.rds"))


## Calculate effect sizes for text: root-to-shoot ####

# Collect MCMC samples for all regression coefficients from root-to-shoot model
ggs(rs_cs_out) %>% 
  filter(substr(Parameter, 1, 4) == "beta" | substr(Parameter, 1, 5) == "mu.al") %>% 
  spread(key = Parameter, value = value) -> rs_betas

# Because weight ic_weight and ln_depth were centered, their mean = 0 so they
# don't need to be included in the predicted equation
rs_betas %>% 
  mutate(frame1.age1.loc1.pred = mu.alpha + `beta[3]`,
         frame1.age1.loc2.pred = mu.alpha + `beta[3]` + `beta[6]`,
         frame1.age2.loc1.pred = mu.alpha + `beta[3]` + `beta[7]`,
         frame1.age2.loc2.pred = mu.alpha + `beta[3]` + `beta[6]` + `beta[7]`,
         frame2.age1.loc1.pred = mu.alpha + `beta[4]`,
         frame2.age1.loc2.pred = mu.alpha + `beta[4]` + `beta[6]`,
         frame2.age2.loc1.pred = mu.alpha + `beta[4]` + `beta[7]`,
         frame2.age2.loc2.pred = mu.alpha + `beta[4]` + `beta[6]` + `beta[7]`,
         frame3.age1.loc1.pred = mu.alpha + `beta[5]`,
         frame3.age1.loc2.pred = mu.alpha + `beta[5]` + `beta[6]`,
         frame3.age2.loc1.pred = mu.alpha + `beta[5]` + `beta[7]`,
         frame3.age2.loc2.pred = mu.alpha + `beta[5]` + `beta[6]` + `beta[7]`,
         frame4.age1.loc1.pred = mu.alpha,
         frame4.age1.loc2.pred = mu.alpha + `beta[6]`,
         frame4.age2.loc1.pred = mu.alpha + `beta[7]`,
         frame4.age2.loc2.pred = mu.alpha + `beta[6]` + `beta[7]`) -> predicted_rs

# Calculate average predicted root-to-shoot for Corn Island provenance
predicted_rs %>% 
  select(contains("loc1")) %>% 
  rowMeans() -> rs_corn

# Calculate average predicted root-to-shoot for Sellman Creek provenance
predicted_rs %>% 
  select(contains("loc2")) %>% 
  rowMeans() -> rs_sellman

# Calculate absolute difference in means
mean(rs_corn - rs_sellman) # 0.1161002
# Calculate 95% credible interval of difference in means
quantile(rs_corn - rs_sellman, c(0.025, 0.975)) # 0.01702895 0.21188719
# Calculate mean percent increase from Sellman to Corn
mean(rs_corn / rs_sellman - 1) # 0.1776576
# Calculate 95% credible interval percent increase from Sellman to Corn
quantile(rs_corn/rs_sellman - 1, c(0.025, 0.975)) # 0.02320939 0.35092684 

# Calculate average predicted root-to-shoot for ancestral cohort
predicted_rs %>% 
  select(contains("age1")) %>% 
  rowMeans() -> rs_ancestral

# Calculate average predicted root-to-shoot for modern cohort
predicted_rs %>% 
  select(contains("age2")) %>% 
  rowMeans() -> rs_modern

# Calculate mean percent decrease from ancestral to modern
mean((rs_ancestral - rs_modern)/rs_ancestral) # 0.08648667
# Calculate 95% quantile for percent decrease from ancestral to modern
quantile((rs_ancestral - rs_modern)/rs_ancestral, c(0.025, 0.975)) # -0.03826322  0.19962476  


## Calculate effect sizes for text: stem width ####

# Pull out MCMC samples from width model
ggs(width_cs_out) %>% 
  filter(substr(Parameter, 1, 4) == "beta" | substr(Parameter, 1, 5) == "mu.al") %>% 
  spread(key = Parameter, value = value) -> width_betas

# Calculate predicted means at each combination of frame, location and age
# (don't need to include continuous covariates because data are centered on 0)
width_betas %>% 
  mutate(frame1.age1.loc1.pred = mu.alpha + `beta[3]`,
         frame1.age1.loc2.pred = mu.alpha + `beta[3]` + `beta[6]`,
         frame1.age2.loc1.pred = mu.alpha + `beta[3]` + `beta[7]`,
         frame1.age2.loc2.pred = mu.alpha + `beta[3]` + `beta[6]` + `beta[7]`,
         frame2.age1.loc1.pred = mu.alpha + `beta[4]`,
         frame2.age1.loc2.pred = mu.alpha + `beta[4]` + `beta[6]`,
         frame2.age2.loc1.pred = mu.alpha + `beta[4]` + `beta[7]`,
         frame2.age2.loc2.pred = mu.alpha + `beta[4]` + `beta[6]` + `beta[7]`,
         frame3.age1.loc1.pred = mu.alpha + `beta[5]`,
         frame3.age1.loc2.pred = mu.alpha + `beta[5]` + `beta[6]`,
         frame3.age2.loc1.pred = mu.alpha + `beta[5]` + `beta[7]`,
         frame3.age2.loc2.pred = mu.alpha + `beta[5]` + `beta[6]` + `beta[7]`,
         frame4.age1.loc1.pred = mu.alpha,
         frame4.age1.loc2.pred = mu.alpha + `beta[6]`,
         frame4.age2.loc1.pred = mu.alpha + `beta[7]`,
         frame4.age2.loc2.pred = mu.alpha + `beta[6]` + `beta[7]`) -> predicted_widths

# Calculate average width for ancestral
predicted_widths %>% 
  select(contains("age1")) %>% 
  rowMeans() -> widths_ancestral

# Calculate average width for modern
predicted_widths %>% 
  select(contains("age2")) %>% 
  rowMeans() -> widths_modern

# Calculate mean percent decrease from ancestral to modern
mean((widths_ancestral - widths_modern)/widths_ancestral) # 0.05896701
# Calculate 95% credible interval of percent decrease from ancestral to modern
quantile((widths_ancestral - widths_modern)/widths_ancestral, c(0.025, 0.975)) # -0.03401524  0.14842040   


## Calculate effect sizes for text: stem height ####

# Pull MCMC samples of regression coefficients from height model
ggs(height_cs_out) %>% 
  filter(substr(Parameter, 1, 4) == "beta" | substr(Parameter, 1, 5) == "mu.al") %>% 
  spread(key = Parameter, value = value) -> height_betas

# Calculate predicted height for all combinations of age, location, and frame.
# Don't need to include continuous covariates because they were centered on 0.
height_betas %>% 
  mutate(frame1.age1.loc1.pred = mu.alpha + `beta[3]`,
         frame1.age1.loc2.pred = mu.alpha + `beta[3]` + `beta[6]`,
         frame1.age2.loc1.pred = mu.alpha + `beta[3]` + `beta[7]`,
         frame1.age2.loc2.pred = mu.alpha + `beta[3]` + `beta[6]` + `beta[7]`,
         frame2.age1.loc1.pred = mu.alpha + `beta[4]`,
         frame2.age1.loc2.pred = mu.alpha + `beta[4]` + `beta[6]`,
         frame2.age2.loc1.pred = mu.alpha + `beta[4]` + `beta[7]`,
         frame2.age2.loc2.pred = mu.alpha + `beta[4]` + `beta[6]` + `beta[7]`,
         frame3.age1.loc1.pred = mu.alpha + `beta[5]`,
         frame3.age1.loc2.pred = mu.alpha + `beta[5]` + `beta[6]`,
         frame3.age2.loc1.pred = mu.alpha + `beta[5]` + `beta[7]`,
         frame3.age2.loc2.pred = mu.alpha + `beta[5]` + `beta[6]` + `beta[7]`,
         frame4.age1.loc1.pred = mu.alpha,
         frame4.age1.loc2.pred = mu.alpha + `beta[6]`,
         frame4.age2.loc1.pred = mu.alpha + `beta[7]`,
         frame4.age2.loc2.pred = mu.alpha + `beta[6]` + `beta[7]`) -> predicted_heights

# Calculate predicted height for Corn Island
predicted_heights %>% 
  dplyr::select(contains("loc1")) %>% 
  rowMeans() -> heights_corn

# Calculate predicted height for Sellman Creek
predicted_heights %>% 
  dplyr::select(contains("loc2")) %>% 
  rowMeans() -> heights_sellman

# Calculate mean percent increase from Corn to Sellman
mean((heights_sellman - heights_corn)/heights_corn) # 0.0297219
# Calculate 95% credible interval percent increase from Corn to Sellman
quantile((heights_sellman - heights_corn)/heights_corn, c(0.025, 0.975)) # -0.01528762  0.07602574   


## Monoculture vs polyculture analysis ####

# Create a dataset with just polycultures
poly_traits <- all_traits %>% filter(diversity == "poly")

# Create function to calculate predicted trait value given additive interactions
# between genotypes
additive_predict <- function(monoculture_ggs){
  # Set up information about polyculture pots (ln_depth and frame). Not actually
  # fitting a model here, just getting information from the model matrix
  mod_poly <- lm(agb ~ ln_depth + frame, data = poly_traits) 
  
  # Get model matrices from linear models (remove intercept)
  M_poly <- as.matrix(model.matrix(mod_poly)[,-1])
  
  # Scale ln_depth (using means and sd from in monoculture models)
  M_poly_sc <- M_poly
  ln_depth_Mono_mean <- mean(mono_traits$ln_depth)
  ln_depth_Mono_sd <- sd(mono_traits$ln_depth)
  M_poly_sc[,"ln_depth"] <- (M_poly[,"ln_depth"] - ln_depth_Mono_mean)/ln_depth_Mono_sd
  
  # Set up storage to hold biomass predictions
  MonoPredict <- matrix(0, nrow = 6666, ncol = 48)
  # Set up storage
  predict <- matrix(0, nrow = 6666, ncol = 4)
  
  for (i in 1:48){
    
    # We need initial weights for each propagule placed in the polyculture pot
    
    # Multiply by 4 to act as if this pot has four identical propagules. We don't
    # want to be predicting with initial weights outside of ranges that we fitted
    # the regression in.
    
    ic_data %>% 
      filter(pot == i) %>% 
      mutate(weight_4 = weight * 4) -> pot_mono
    
    # Pull covariate information
    M_poly_sc[i,] -> pot_poly
    
    # Pull out all random intercept (genotype) random intercepts
    monoculture_ggs %>%
      filter(substr(Parameter, 1, 5) == "alpha") %>% 
      spread(key = Parameter, value = value) %>% 
      dplyr::select(contains("alpha")) -> alphas_ggs
    
    alphas_df <- as.data.frame(alphas_ggs)
    # This puts genotypes in alphabetical order which matches up with how they
    # were fit in monoculture models. This checks out with raw data too
    colnames(alphas_df) <- levels(factor(mono_traits$genotype))
    
    # Get regression coefficient for ln_depth from monoculture model
    monoculture_ggs %>% 
      filter(Parameter == "beta[1]") %>% 
      pull(value) -> beta_lndepth
    
    # Get regression coefficient for ic_weight
    monoculture_ggs %>% 
      filter(Parameter == "beta[2]") %>% 
      pull(value) -> beta_icweight
    
    # Get first frame regression coefficient
    monoculture_ggs %>% 
      filter(Parameter == "beta[3]") %>% 
      pull(value) -> beta_frame1
    # Get second frame regression coefficient
    monoculture_ggs %>% 
      filter(Parameter == "beta[4]") %>% 
      pull(value) -> beta_frame2
    # Get last frame regression coefficient
    monoculture_ggs %>% 
      filter(Parameter == "beta[5]") %>% 
      pull(value) -> beta_frame3
    
    # Create regression with frame, ic_weight, and depth info from the
    # polyculture and use regression coefficients from the monoculture model
    for (j in 1:4){
      genotype_name <- pot_mono[j, "genotype"]
      # Scale initial condition weight using scaling parameters from monoculture models
      ic_weight_scaled <- (pot_mono[j, "weight_4"] - mean(mono_traits$ic_weight)) / sd(mono_traits$ic_weight)
      predict[,j] <- alphas_df[, genotype_name] + 
        beta_icweight * ic_weight_scaled + beta_lndepth * pot_poly[1] +
        beta_frame1 * pot_poly[2] + beta_frame2 * pot_poly[3] + beta_frame3 * pot_poly[4]
    }
    
    # Then take the row means to get the predicted biomass if this was all additive
    MonoPredict[,i] <- rowMeans(predict, na.rm = T)
  }
  return(list(MonoPredict = MonoPredict))
}

# Run additive function for each trait
agb_additive <- additive_predict(ggs(agb_out))
bgb_additive <- additive_predict(ggs(bgb_out))
height_additive <- additive_predict(ggs(height_out))
width_additive <- additive_predict(ggs(width_out))
rs_additive <- additive_predict(ggs(rs_out))
density_additive <- additive_predict(ggs(density_out))
beta_additive <- additive_predict(ggs(beta_out))

# Save additive output for later plots
saveRDS(agb_additive, here("outputs/monoculture_polyculture/","monopoly_agb_add.rds"))
saveRDS(bgb_additive, here("outputs/monoculture_polyculture/","monopoly_bgb_add.rds"))
saveRDS(height_additive, here("outputs/monoculture_polyculture/","monopoly_height_add.rds"))
saveRDS(width_additive, here("outputs/monoculture_polyculture/","monopoly_width_add.rds"))
saveRDS(rs_additive, here("outputs/monoculture_polyculture/","monopoly_rs_add.rds"))
saveRDS(density_additive, here("outputs/monoculture_polyculture/","monopoly_density_add.rds"))
saveRDS(beta_additive, here("outputs/monoculture_polyculture/","monopoly_beta_add.rds"))

## Calculate differences between observed trait values and additive predictions
## for each pot (j) at each iteration (i)

# Create storage to hold differences between predicted and observed traits. This
# should be the same dimensions as the additive predictions matrix
pred_obs_diff <- matrix(NA, nrow = nrow(biomass_additive$MonoPredict),
                        ncol = ncol(biomass_additive$MonoPredict))

# Loop through polyculture pots (n = 48) to get differences at each iteration.
# This is for figure 2a.
average_difference <- function(trait, additive_samples){
  for (j in 1:48){
    pot_trait_temp <- poly_traits %>%
      filter(pot == j) %>%
      select(trait) %>% as.numeric()
    if(trait == "beta" & is.na(pot_trait_temp)){
      pot_trait_temp <- 10000
    }
    # Observed - Predicted (positive values mean that what we observed is greater
    # than the additive expectation, so POSITIVE = FACTILITATION; and negative
    # values mean that what we observed was less than the additive expectation, so
    # NEGATIVE = COMPETITION)
    if(trait == "density"){
      pred_obs_diff[,j] <- log(pot_trait_temp) - additive_samples$MonoPredict[,j]
    }else{
      pred_obs_diff[,j] <- pot_trait_temp - additive_samples$MonoPredict[,j]
    }
     
  }
  # Then take a mean across all 48 pots for each iteration to get the mean
  # effect (for beta make sure to skip pot 39)
  if(trait == "beta"){
    average_difference <- rowMeans(pred_obs_diff[,c(1:38,40:48)])
  }else{
    average_difference <- rowMeans(pred_obs_diff)
  }
  return(average_difference)
}

## Differences by age cohort: Calculate differences at each iteration for each
## pot. This is for figure 2b and supplemental figs.
mean_difference_bypot <- function(trait, additive_samples){
  # Predicted values for each polyculture pot (row = iterations, col = pots)
  observed <- pull(poly_traits[,trait])
  
  # Create matrix to hold differences between observed and predicted for each pot
  # at each iteration
  difference <- matrix(0, nrow = 6666, ncol = 48)
  
  if(trait == "density"){
    for (i in 1:6666){
      difference[i,] <- log(observed) - additive_samples$MonoPredict[i,]
    }
  }else{
    for (i in 1:6666){
      difference[i,] <- observed - additive_samples$MonoPredict[i,]
    }
  }
  
  # Get average difference across pots for each iteration
  avg_difference <- data.frame(x =colMeans(difference, na.rm = T))
  
  return(avg_difference)
}

# Calculate average difference for each pot for each trait
diffs_biomass <- mean_difference_bypot("agb", biomass_additive)
diffs_density <- mean_difference_bypot("density", density_additive)
diffs_height<- mean_difference_bypot("mean_tot_height", height_additive)
diffs_width <- mean_difference_bypot("mean_mid_width", width_additive)
diffs_bgb <- mean_difference_bypot("bgb", bgb_additive)
diffs_rs <- mean_difference_bypot("rs", rs_additive)
diffs_beta <- mean_difference_bypot("beta", beta_additive)

# Create a data frame with all average differences and age cohort information
tibble(poly_traits) %>%
  arrange(pot) %>%
  mutate(`aboveground biomass (g)` = diffs_biomass$x,
         `stem density` = diffs_density$x,
         `mean stem height (cm)` = diffs_height$x,
         `mean stem width (mm)` = diffs_width$x,
         `belowground biomass (g)` = diffs_bgb$x,
         `root:shoot ratio` = diffs_rs$x,
         `root distribution parameter` = diffs_beta$x) %>%
  mutate(age = case_when(pot %in% 1:18 ~ "ancestral",
                         pot %in% 19:30 ~ "modern",
                         T ~ "mix")) %>%
  dplyr::select(age, `aboveground biomass (g)`, `stem density`, `mean stem height (cm)`, `mean stem width (mm)`,
                `belowground biomass (g)`, `root:shoot ratio`, `root distribution parameter`) -> diffs_by_age

# Save table for later
saveRDS(diffs_by_age, here("outputs/monoculture_polyculture/", "diffs_by_age.rds"))

# Check to see if there are significant differences by age group
abg_mod <- lm(`aboveground biomass (g)` ~ age, data = diffs_by_age)
anova(abg_mod) # ns

bgb_mod <- lm(`belowground biomass (g)` ~ age, data = diffs_by_age)
anova(bgb_mod) # ns

density_mod <- lm(`stem density` ~ age, data = diffs_by_age)
anova(density_mod) # ns

beta_mod <- lm(`root distribution parameter` ~ age, data = diffs_by_age)
anova(beta_mod) # .
emmeans::emmeans(beta_mod, ~ age) # mix has lower beta than would be expected based on additive

rs_mod <- lm(`root:shoot ratio` ~ age, data = diffs_by_age)
anova(rs_mod) # ns

height_mod <- lm(`mean stem height (cm)` ~ age, data = diffs_by_age)
anova(height_mod) # ns

width_mod <- lm(`mean stem width (mm)` ~ age, data = diffs_by_age)
anova(width_mod) # ns


## Cohort Marsh Equilibrium Model Simulations

library(rCTM); library(tidyverse); library(ggmcmc)
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

## Preliminaries ####

# Parameter estimates for bMax and root:shoot from Blue Genes 2019 experiment
source(here("chp1/code/MS_code/mem_params_BlueGenes.R"))

# Create function to collect mean and sd from each bayesian regression model
get_cv <- function(coda_object){
  ggs(coda_object) %>% 
    filter(Parameter %in% c("sigma.int", "mu.alpha")) %>% 
    group_by(Parameter) %>% 
    summarize(mean = mean(value)) %>% 
    spread(key = Parameter, value = mean) %>% 
    mutate(cv = sigma.int / mu.alpha) -> out
  
  return(out)
}

# Collect mean and sd of each of the traits: aboveground biomass, root:shoot
# ratio, and root distribution parameter
biomass_sum <- get_cv(readRDS("chp1/results/biomass_monomodel.rds"))
rs_sum <- get_cv(readRDS("chp1/results/rs_monomodel.rds"))
beta_sum <- get_cv(readRDS("chp1/results/beta_monomodel.rds"))

mean_vec <- c(BlueGenes_Params_forsim$bMax, BlueGenes_Params_forsim$rootToShoot, beta_sum$mu.alpha)
#mean_vec <- c(0.1, 2, beta_sum$mu.alpha)

# Scale biomass to g/cm^2
# Pots are 7.62^2*pi cm2
pot_area <- 7.62^2*pi
biomass_sum %>% 
  mutate(mu.alpha.scaled = mu.alpha / pot_area,
         sigma.int.scaled = sigma.int / pot_area) -> biomass_sum

## All variation - vary due to genotype ####
# Create covariance matrix if we assume that variance is the same, but means are
# different (conservative assumption)
mono_all <- readRDS("chp1/results/mono_all_derived.rds")

# Maximum biomass estimates from blue genes
biomass_scale <- mean_vec[1] / biomass_sum$mu.alpha.scaled
rs_scale <- mean_vec[2] / rs_sum$mu.alpha

var_biomass <- (biomass_sum$sigma.int.scaled * biomass_scale)^2
var_rs <- (rs_sum$sigma.int * rs_scale)^2
var_beta <- beta_sum$sigma.int^2 # keep beta the same because it's not off of what we might expect


# var_biomass <- biomass_sum$sigma.int.scaled^2
# var_rs <- rs_sum$sigma.int^2
# var_beta <- beta_sum$sigma.int^2

covar_biomass_rs <- biomass_sum$sigma.int.scaled * rs_sum$sigma.int * cor(mono_all$total_biomass, mono_all$rs)
covar_biomass_beta <- biomass_sum$sigma.int.scaled * beta_sum$sigma.int * cor(mono_all$total_biomass, mono_all$beta)
covar_rs_beta <- rs_sum$sigma.int * beta_sum$sigma.int * cor(mono_all$rs, mono_all$beta)

covar_matrix <- matrix(c(var_biomass, covar_biomass_rs, covar_biomass_beta,
                         covar_biomass_rs, var_rs, covar_rs_beta,
                         covar_biomass_beta, covar_rs_beta, var_beta),
                       nrow = 3, byrow = T)

# Take 1000 draws from multivariate normal distribution
samples1 <- mvtnorm::rmvnorm(1000, mean = mean_vec, covar_matrix)

# Translate beta values to maximum rooting depth
depth_interval <- seq(0,50,length.out = 1000)
rooting_depth <- NULL

for(i in 1:nrow(samples1)){
  temp_beta <- samples1[i,3]
  temp_cumulative <- 1 - temp_beta ^ depth_interval
  depth_95 <- which.min(abs(temp_cumulative - 0.95))
  rooting_depth[i] <- depth_interval[depth_95]
}

# Put all data for simulations into a data frame
tibble(`aboveground biomass (g)` = samples1[,1],
       `root:shoot ratio` = samples1[,2],
       `maximum rooting depth (cm)` = rooting_depth) -> for_MEM1

saveRDS(for_MEM1,"chp1/results/randomdrawsMEM.rds")

# Create storage for MEM model runs
n_runs <- 100
run_store1 <- matrix(NA, nrow = n_runs, ncol = 100)
carbon_store <- matrix(NA, nrow = n_runs, ncol = 100)
flood_store <- matrix(NA, nrow = n_runs, ncol = 100)

for (i in 1:n_runs){
  mem_out <- runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = for_MEM1$`aboveground biomass (g)`[i], 
                          zVegMin=BlueGenes_Params_forsim$zVegMin, zVegMax=BlueGenes_Params_forsim$zVegMax, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = for_MEM1$`root:shoot ratio`[i],
                          rootTurnover=0.5, rootDepthMax=for_MEM1$`maximum rooting depth (cm)`[i], omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store1[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  flood_store[i,] <- mem_out$annualTimeSteps$meanHighWater - mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store[i,]
  print(i)
}

par(mar = c(4,4,2,2))
plot(run_store1[1,1:80], type = "n", ylim = c(22.6,35), xlab = "time forward (years)",
     ylab = "marsh elevation (cm)")
for(i in 1:n_runs){
  lines(run_store1[i,1:80], col = rgb(0,0,0,0.2))
}

# Calculate CIs around average accretion rate
init_elev <- 22.6
avg_accretion_rates1 <- (run_store1[,80] - init_elev) / 80
avg_acc_rate_ci1 <- quantile(avg_accretion_rates1, c(0.025, 0.5, 0.975))
avg_acc_rate_ci1[3]/avg_acc_rate_ci1[1] #109% increase (more than double)

# Calculate CIs around carbon accumulation rate
avg_C_accum_rate <- (carbon_store[,80] - carbon_store[,1]) / 80
avg_C_accum_rate_ci <- quantile(avg_C_accum_rate, c(0.025, 0.5, 0.975))
avg_C_accum_rate_ci[3]/avg_C_accum_rate_ci[1] # 158% increase (more than double)

# Convert C accumulation rates to metric tons per hectare (x metric tons / hectare conversion factors)
avg_C_accum_rate * 1e-6 / 1e-8

## All variation - vary by cohort mean ####
cs_all <- mono_all %>% filter(location %in% c('corn', "sellman"))

# Root:shoot
rs_mod_formeans <- lmer(rs ~ frame + ic_weight + ln_depth + location + age + (1|genotype), data = cs_all)
means_rs <- summary(emmeans(rs_mod_formeans, ~ location:age))$emmean
bg_rs_scale <- BlueGenes_Params_forsim$rootToShoot / mean(means_rs)
root_shoot_cohort_forMEM <- means_rs * bg_rs_scale

# Biomass
biomass_mod_formeans <- lmer(total_biomass ~ ic_weight + frame + ln_depth + location + age + (1|genotype), data = cs_all)
means_biomass <- summary(emmeans(biomass_mod_formeans, ~ location:age))$emmean / pot_area
bg_biomass_scale <- BlueGenes_Params_forsim$bMax / mean(means_biomass)
biomass_cohort_forMEM <- means_biomass * bg_biomass_scale

# Root distribution parameter
beta_mod_formeans <- lm(beta ~ frame + ic_weight + ln_depth + location + age, data = cs_all)
means_betas <- summary(emmeans(beta_mod_formeans, ~ location:age))$emmean

# Translate beta values to maximum rooting depth
depth_interval <- seq(0,50,length.out = 1000)
rooting_depth <- NULL

for(i in 1:length(means_betas)){
  temp_beta <- means_betas[i]
  temp_cumulative <- 1 - temp_beta ^ depth_interval
  depth_95 <- which.min(abs(temp_cumulative - 0.95))
  rooting_depth[i] <- depth_interval[depth_95]
}

betas_cohort_forMEM <- rooting_depth

# Run model for each cohort
run_store_cohort <- matrix(NA, nrow = 4, ncol = 100)
carbon_store_cohort <- matrix(NA, nrow = 4, ncol = 100)
flood_store_cohort <- matrix(NA, nrow = 4, ncol = 100)

for (i in 1:4){
  mem_out <- runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = biomass_cohort_forMEM[i], 
                          zVegMin=BlueGenes_Params_forsim$zVegMin, zVegMax=BlueGenes_Params_forsim$zVegMax, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = root_shoot_cohort_forMEM[i],
                          rootTurnover=0.5, rootDepthMax=rooting_depth[i], omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_cohort[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  flood_store_cohort[i,] <- mem_out$annualTimeSteps$meanHighWater - mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_cohort[i,]
  print(i)
}

# Calculate CIs around average accretion rate
init_elev <- 22.6
avg_accretion_rates_cohort <- (run_store_cohort[,80] - init_elev) / 80
avg_accretion_rates_cohort[1]/avg_accretion_rates_cohort[4] # 24% increase in accretion rate from lowest to highest

# Calculate CIs around carbon accumulation rate
avg_C_accum_rate_cohort <- (carbon_store_cohort[,80] - carbon_store_cohort[,1]) / 80
avg_C_accum_rate_cohort[1]/avg_C_accum_rate_cohort[4] # 32% increase in carbon accumulation rate from lowest to highest

# # Convert C accumulation rates to metric tons per hectare (x metric tons / hectare conversion factors)
# avg_C_accum_rate_cohort * 1e-6 / 1e-8

tibble(location = c("corn", "sellman", "corn", "sellman"),
       age = c("ancestral", "ancestral", "modern", "modern"),
       acc_v = avg_accretion_rates_cohort * 10,
       acc_C = avg_C_accum_rate_cohort * 1e-6 / 1e-8) -> cohort_summary

elevation_store_df <- as.data.frame(run_store1)
elevation_store_withcohorts <- rbind(elevation_store_df[,1:80], run_store_cohort[,1:80])
colnames(elevation_store_withcohorts) <- 1:80

tibble(elevation_store_withcohorts) %>% 
  mutate(iteration = 1:104) %>% 
  gather(key = year, value = value, `1`:`80`) %>% 
  mutate(year = as.numeric(year) + 2020) %>% 
  mutate(color_code = case_when(iteration == 101 ~ 1,
                                iteration == 102 ~ 2,
                                iteration == 103 ~ 3,
                                iteration == 104 ~ 4,
                                T ~ 0),
         size_code = case_when(iteration > 100 ~ "big",
                               T ~ "small")) %>% 
  ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) + 
  geom_line(aes(alpha = size_code)) +
  ylab("marsh elevation (cm NAVD88)") +
  scale_color_manual(values = c("gray11", colors[1], colors[4], colors[1], colors[4])) +
  scale_size_manual(values = c(1.5,0.8)) +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  geom_point(aes(x = 2100, y = run_store_cohort[1,80]), color = colors[1], size = 3) +
  geom_point(aes(x = 2100, y = run_store_cohort[2,80]), color = colors[4], size = 3) +
  geom_point(aes(x = 2100, y = run_store_cohort[3,80]), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = 2100, y = run_store_cohort[4,80]), color = colors[4], size = 3, shape = 17) +
  theme(legend.position = "none") +
  ylim(22,33.5)-> Fig4_panelA

tibble(acc_rate = avg_accretion_rates1*10) %>% 
  ggplot(aes(x = acc_rate)) +
  geom_histogram(bins=10, color = "black", fill = "white") +
  xlab(expression(paste("vertical accretion rate (mm ",yr^-1,")"))) +
  geom_point(aes(x = cohort_summary$acc_v[1], y = 0), color = colors[1], size = 3) +
  geom_point(aes(x = cohort_summary$acc_v[2], y = 0), color = colors[4], size = 3) +
  geom_point(aes(x = cohort_summary$acc_v[3], y = 0), color = colors[1], size = 3, shape = 17)+
  geom_point(aes(x = cohort_summary$acc_v[4], y = 0), color = colors[4], size = 3, shape = 17) -> Fig4_panelB

tibble(acc_rate = avg_C_accum_rate * 1e-6 / 1e-8) %>% 
  ggplot(aes(x = acc_rate)) +
  geom_histogram(bins = 10, fill = "white", color = "black") +
  xlab(expression(paste("carbon accumulation rate (t C ", ha^-1, yr^-1,")"))) +
  geom_point(aes(x = cohort_summary$acc_C[1], y = 0), color = colors[1], size = 3) +
  geom_point(aes(x = cohort_summary$acc_C[2], y = 0), color = colors[4], size = 3) +
  geom_point(aes(x = cohort_summary$acc_C[3], y = 0), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = cohort_summary$acc_C[4], y = 0), color = colors[4], size = 3, shape = 17) -> Fig4_panelC 

Fig4_panelsBC <- cowplot::plot_grid(Fig4_panelB, Fig4_panelC, nrow = 2, labels = c("b", "c"))
Fig4_panelA_label <- cowplot::plot_grid(Fig4_panelA, labels = "a")

##
# Calculations for in-text
##

# Differences in final marsh elevation due to genotype
mean(run_store1[,80]) # 29.13397
quantile(run_store1[,80], c(0.025, 0.975)) # 27.34588 32.22798 

# Vertical accretion rates
init_elev <- 22.6
avg_accretion_rates1 <- (run_store1[,80] - init_elev) / 80
mean(avg_accretion_rates1)
quantile(avg_accretion_rates1, c(0.025, 0.975))

# Carbon accumulation rates
avg_C_accum_rate <- (carbon_store[,80] - carbon_store[,1]) / 80
mean(avg_C_accum_rate)
quantile(avg_C_accum_rate, c(0.025, 0.975))
quantile(avg_C_accum_rate, c(0.975)) - quantile(avg_C_accum_rate, c(0.025))

# % diff accretion rate for ecotypes and age cohorts
# Increase accretion rate from sellman to corn
mean(cohort_summary$acc_v[c(1,3)]) / mean(cohort_summary$acc_v[c(2,4)])
# Decrease in accretion rate from ancestral to modern
(mean(cohort_summary$acc_v[c(1,2)]) - mean(cohort_summary$acc_v[c(3,4)])) / mean(cohort_summary$acc_v[c(1,2)])
# Increase in accretion rate from modern to ancestral
mean(cohort_summary$acc_v[c(1,2)]) / mean(cohort_summary$acc_v[c(3,4)])

# Same for C accumulation
# Increase accretion rate from sellman to corn
mean(cohort_summary$acc_C[c(1,3)]) / mean(cohort_summary$acc_C[c(2,4)])
# Decrease in accretion rate from ancestral to modern
(mean(cohort_summary$acc_C[c(1,2)]) - mean(cohort_summary$acc_C[c(3,4)])) / mean(cohort_summary$acc_C[c(1,2)])
# Increase in accretion rate from modern to ancestral
mean(cohort_summary$acc_C[c(1,2)]) / mean(cohort_summary$acc_C[c(3,4)])


## Allow just agb to vary - vary due to genotype ####

# Take 100 draws from normal distribution just for biomass
samples_Biomass <-rnorm(100, mean = mean_vec[1], sqrt(var_biomass))

# Translate beta values to maximum rooting depth
depth_interval <- seq(0,50,length.out = 1000)
rooting_depth <- NULL

temp_beta <- mean_vec[3]
temp_cumulative <- 1 - temp_beta ^ depth_interval
depth_95 <- which.min(abs(temp_cumulative - 0.95))
rooting_depth <- depth_interval[depth_95]

# Put all data for simulations into a data frame
tibble(`aboveground biomass (g)` = samples_Biomass,
       `root:shoot ratio` = mean_vec[2],
       `maximum rooting depth (cm)` = rooting_depth) -> for_MEM2

saveRDS(for_MEM2,"chp1/results/randomdrawsMEM_biomass_var.rds")

# Create storage for MEM model runs
n_runs <- 100
run_store2 <- matrix(NA, nrow = n_runs, ncol = 100)
carbon_store2 <- matrix(NA, nrow = n_runs, ncol = 100)
flood_store2 <- matrix(NA, nrow = n_runs, ncol = 100)

for (i in 1:n_runs){
  mem_out <- runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = for_MEM2$`aboveground biomass (g)`[i], 
                          zVegMin=BlueGenes_Params_forsim$zVegMin, zVegMax=BlueGenes_Params_forsim$zVegMax, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = for_MEM2$`root:shoot ratio`[1],
                          rootTurnover=0.5, rootDepthMax=for_MEM2$`maximum rooting depth (cm)`[1], omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store2[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  flood_store2[i,] <- mem_out$annualTimeSteps$meanHighWater - mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store2[i,]
  print(i)
}

par(mar = c(4,4,2,2))
plot(run_store2[1,1:80], type = "n", ylim = c(22.6,35), xlab = "time forward (years)",
     ylab = "marsh elevation (cm)")
for(i in 1:length(samples_Biomass)){
  lines(run_store2[i,1:80], col = rgb(0,0,0,0.1))
}

# Calculate CIs around average accretion rate
init_elev <- 22.6
avg_accretion_rates2 <- (run_store2[,80] - init_elev) / 80
avg_acc_rate_ci2 <- quantile(avg_accretion_rates2, c(0.025, 0.5, 0.975))
avg_acc_rate_ci2[3]/avg_acc_rate_ci2[1] # 56% increase

# Calculate CIs around carbon accumulation rate
avg_C_accum_rate2 <- (carbon_store2[,80] - carbon_store2[,1]) / 80
avg_C_accum_rate_ci2 <- quantile(avg_C_accum_rate2, c(0.025, 0.5, 0.975))
avg_C_accum_rate_ci2[3]/avg_C_accum_rate_ci2[1] # 76.9% increase 

## Allow just agb to vary - vary due to cohort ####
run_store_cohort2 <- matrix(NA, nrow = 4, ncol = 100)
carbon_store_cohort2 <- matrix(NA, nrow = 4, ncol = 100)
flood_store_cohort2 <- matrix(NA, nrow = 4, ncol = 100)

# For root:shoot and rooting depth, take means across cohorts
for (i in 1:4){
  mem_out <- runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = biomass_cohort_forMEM[i], 
                          zVegMin=BlueGenes_Params_forsim$zVegMin, zVegMax=BlueGenes_Params_forsim$zVegMax, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                          rootTurnover=0.5, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_cohort2[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  flood_store_cohort2[i,] <- mem_out$annualTimeSteps$meanHighWater - mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_cohort2[i,]
  print(i)
}

# Quick plot
par(mar = c(4,4,2,2))
plot(run_store_cohort2[1,1:80], type = "n", ylim = c(22.6,30), xlab = "time forward (years)",
     ylab = "marsh elevation (cm)")
for(i in 1:4){
  lines(run_store_cohort2[i,1:80], col = c(colors[1], colors[4], colors[1], colors[4])[i])
}

# Calculate CIs around average accretion rate
init_elev <- 22.6
avg_accretion_rates_cohort2 <- (run_store_cohort2[,80] - init_elev) / 80
avg_accretion_rates_cohort2[2]/avg_accretion_rates_cohort2[3] # 19.7% increase in accretion rate from lowest to highest

# Calculate CIs around carbon accumulation rate
avg_C_accum_rate_cohort2 <- (carbon_store_cohort2[,80] - carbon_store_cohort2[,1]) / 80
avg_C_accum_rate_cohort2[2]/avg_C_accum_rate_cohort2[3] # 26.0% increase in carbon accumulation rate from lowest to highest

tibble(location = c("corn", "sellman", "corn", "sellman"),
       age = c("ancestral", "ancestral", "modern", "modern"),
       acc_v = avg_accretion_rates_cohort2 * 10,
       acc_C = avg_C_accum_rate_cohort2 * 1e-6 / 1e-8) -> cohort_summary

elevation_store_df2 <- as.data.frame(run_store2)
elevation_store_withcohorts2 <- rbind(elevation_store_df2[,1:80], run_store_cohort2[,1:80])
colnames(elevation_store_withcohorts2) <- 1:80

tibble(elevation_store_withcohorts2) %>% 
  mutate(iteration = 1:104) %>% 
  gather(key = year, value = value, `1`:`80`) %>% 
  mutate(year = as.numeric(year) + 2020) %>% 
  mutate(color_code = case_when(iteration == 101 ~ 1,
                                iteration == 102 ~ 2,
                                iteration == 103 ~ 3,
                                iteration == 104 ~ 4,
                                T ~ 0),
         size_code = case_when(iteration > 100 ~ "big",
                               T ~ "small")) %>% 
  ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) + 
  geom_line(aes(alpha = size_code)) +
  ylab("marsh elevation (cm NAVD88)") +
  scale_color_manual(values = c("gray67", colors[1], colors[4], colors[1], colors[4])) +
  scale_size_manual(values = c(1.5,0.8)) +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  geom_point(aes(x = 2100, y = run_store_cohort2[1,80]), color = colors[1], size = 3) +
  geom_point(aes(x = 2100, y = run_store_cohort2[2,80]), color = colors[4], size = 3,) +
  geom_point(aes(x = 2100, y = run_store_cohort2[3,80]), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = 2100, y = run_store_cohort2[4,80]), color = colors[4], size = 3, shape = 17) +
  theme(legend.position = "none") +
  ylim(22,33.5)-> Fig4_panelD

## Last panel -- compare variation for agb + bgb and abg only models ####

quantile_agb_only <- quantile(run_store2[,80], c(0.025, 0.975))
quantile_agb_bgb <- quantile(run_store1[,80], c(0.025, 0.975))

tibble(y = c(run_store2[,80], run_store1[,80]),
       x = rep(c("agb only", "agb + bgb"), each = length(run_store1[,80]))) %>% 
  ggplot(aes(x = x, y = y, fill = x)) +
  geom_violin(draw_quantiles = c(0.025,0.5, 0.975), size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("gray11", "gray67")) +
  theme(legend.position = "none") +
  ylab("elevation at t = 80 (cm NAVD88)") +
  xlab("scenario") -> Fig4_panelE

## Bring all plots together ####
Fig4_panelsDE <- cowplot::plot_grid(Fig4_panelD, Fig4_panelE, nrow = 1,
                                    labels = c("d", "e"), rel_widths = c(3,2))
Fig4_panelsABC <- cowplot::plot_grid(Fig4_panelA_label, Fig4_panelsBC, rel_widths = c(3,2))

Fig4_panelsABCDE <- cowplot::plot_grid(Fig4_panelsABC, Fig4_panelsDE, nrow = 2)

png("chp1/plots/MS_plots/Vahsen_Fig4.png", height = 6.8, width = 8, units = "in", res = 300)
Fig4_panelsABCDE
dev.off()
