# Vahsen et al. Belowground Evolution main script

# This script runs the main analyses and sources additional functions within the
# main_code folder

## Uncomment lines below to install rCMEM (Cohort Marsh Equilibrium Model) ####

# Uninstall and reinstall developer branch from GitHub

# # 1. If rCMEM is loaded and in the memory, forget rCMEM
# if ("rCMEM" %in% (.packages())){
#   detach("package:rCMEM", unload=TRUE)
# }
# 
# # 2. If remotes is not already installed, install it
# if (! ("remotes" %in% installed.packages())) {
#   install.packages("remotes")
# }
# 
# # 3. Install package from developer branch of GitHub
# devtools::install_github("https://github.com/tilbud/rCMEM/tree/JimH-dev")

## Load libraries ####
library(tidyverse); library(ggmcmc); library(rjags); library(lme4)
library(here); library(rCMEM); library(mvtnorm); library(emmeans)

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

# Set seed - this is optional. Setting seed to 1234 will reproduce exact results
# as reported in the manuscript
seed <- 1234

## Fit all basic trait models ####

## aboveground biomass
agb_out <- run_models(mono_traits, "agb", model_template_mono, diag_plot = F, seed = seed)
## stem density
density_out <- run_models(mono_traits, "density", model_template_mono, diag_plot = F, seed = seed)
# Fails to set monitor for sigma.res because there is no residual variation in
# Poisson regression (aka this warning is ok and expected)
## stem height
height_out <- run_models(mono_traits, "mean_tot_height", model_template_mono, diag_plot = F, seed = seed)
## stem width
width_out <- run_models(mono_traits, "mean_mid_width", model_template_mono, diag_plot = F, seed = seed)
## belowground biomass
bgb_out <- run_models(mono_traits, "bgb", model_template_mono, diag_plot = F, seed = seed)
## root-to-shoot ratio
rs_out <- run_models(mono_traits, "rs", model_template_mono, diag_plot = F, seed = seed)
## beta (belowground biomass distribution parameter)
beta_out <- run_models(mono_traits, "beta", model_template_mono, diag_plot = F, seed = seed)

# Save as R data file objects to use later for plotting and subsequent analyses
saveRDS(agb_out, here("outputs/monoculture_models", "agb_monomodel.rds"))
saveRDS(rs_out, here("outputs/monoculture_models", "rs_monomodel.rds"))
saveRDS(bgb_out, here("outputs/monoculture_models", "bgb_monomodel.rds"))
saveRDS(width_out, here("outputs/monoculture_models", "width_monomodel.rds"))
saveRDS(height_out, here("outputs/monoculture_models", "height_monomodel.rds"))
saveRDS(density_out, here("outputs/monoculture_models", "density_monomodel.rds"))
saveRDS(beta_out, here("outputs/monoculture_models", "beta_monomodel.rds"))

# Calculate ICC for all traits
calculate_icc <- function(coda_object){
  ggs(coda_object) %>% 
    filter(Parameter %in% c("sigma.res", "sigma.int")) %>% 
    spread(key = Parameter, value = value) %>% 
    mutate(ICC = sigma.int^2 / (sigma.int^2 + sigma.res^2)) %>% 
    summarize(mean = mean(ICC),
              lower = quantile(ICC, 0.025),
              upper = quantile(ICC, 0.975))
}

calculate_icc(agb_out)
# mean    lower upper
# 0.230 0.00211 0.578
calculate_icc(width_out)
# mean  lower upper
# 0.440 0.0484 0.760
calculate_icc(height_out)
# mean     lower upper
# 0.159 0.00140 0.499
calculate_icc(rs_out)
# mean lower upper
# 0.691 0.437 0.877
calculate_icc(bgb_out)
# mean lower upper
# 0.516 0.157 0.788
calculate_icc(beta_out)
# mean lower upper
# 0.495 0.157 0.773

# Use Nakagawa et al. equation to calculate ICC for poisson regression
ggs(density_out) %>% 
  spread(key = Parameter, value = value) %>%
  mutate(sigma2.res = log(1 + 1/mean(mono_traits$density))) %>% 
  mutate(sigma2.int = sigma.int^2) %>% 
  mutate(ICC = sigma2.int / (sigma2.int + sigma2.res)) %>% 
  summarize(mean = mean(ICC),
            lower = quantile(ICC, 0.025),
            upper = quantile(ICC, 0.975))
# mean  lower upper
# 0.328 0.0335 0.643

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
agb_cs_out <- run_models(cs_traits, "agb", model_template_cs, diag_plot = F, seed = seed)
## stem density
density_cs_out <- run_models(cs_traits, "density", model_template_cs, diag_plot = F, seed = seed)
# Fails to set monitor for sigma.res because there is no residual variation in
# Poisson regression (aka this warning is ok and expected)
## stem height
height_cs_out <- run_models(cs_traits, "mean_tot_height", model_template_cs, diag_plot = F, seed = seed)
## stem width
width_cs_out <- run_models(cs_traits, "mean_mid_width", model_template_cs, diag_plot = F, seed = seed)
## belowground biomass
bgb_cs_out <- run_models(cs_traits, "bgb", model_template_cs, diag_plot = F, seed = seed)
## root-to-shoot ratio
rs_cs_out <- run_models(cs_traits, "rs", model_template_cs, diag_plot = F, seed = seed)
## beta (belowground biomass distribution parameter)
beta_cs_out <- run_models_fixed(cs_traits, "beta", model_template_cs, diag_plot = F, seed = seed)
# Fails to set trace monitor for sigma.int and mu.alpha because this model does
# not have a random intercept for genotype (and is thus just a fixed effects
# model). Because provenance and age cohort explain so much of the variation,
# between-genotype variation is estimated to be zero and was thus dropped from
# the model.

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
mean(rs_corn - rs_sellman) # 0.1130211
# Calculate 95% credible interval of difference in means
quantile(rs_corn - rs_sellman, c(0.025, 0.975)) # 0.01418939 0.20940235 
# Calculate mean percent increase from Sellman to Corn
mean(rs_corn / rs_sellman - 1) # 0.1723259
# Calculate 95% credible interval percent increase from Sellman to Corn
quantile(rs_corn/rs_sellman - 1, c(0.025, 0.975)) # 0.0191938 0.3444169 

# Calculate average predicted root-to-shoot for ancestral cohort
predicted_rs %>% 
  select(contains("age1")) %>% 
  rowMeans() -> rs_ancestral

# Calculate average predicted root-to-shoot for modern cohort
predicted_rs %>% 
  select(contains("age2")) %>% 
  rowMeans() -> rs_modern

# Calculate mean percent decrease from ancestral to modern
mean((rs_ancestral - rs_modern)/rs_ancestral) # 0.08258125
# Calculate 95% quantile for percent decrease from ancestral to modern
quantile((rs_ancestral - rs_modern)/rs_ancestral, c(0.025, 0.975)) # -0.05186102  0.19569098  


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
mean((widths_ancestral - widths_modern)/widths_ancestral) # 0.05639298
# Calculate 95% credible interval of percent decrease from ancestral to modern
quantile((widths_ancestral - widths_modern)/widths_ancestral, c(0.025, 0.975)) # -0.03271603  0.14385558 

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
mean((heights_sellman - heights_corn)/heights_corn) # 0.03032149
# Calculate 95% credible interval percent increase from Corn to Sellman
quantile((heights_sellman - heights_corn)/heights_corn, c(0.025, 0.975)) # -0.01503261  0.07468390   

# Repeat the same process for age
# Calculate predicted height for ancestral
predicted_heights %>% 
  dplyr::select(contains("age1")) %>% 
  rowMeans() -> heights_anc

# Calculate predicted height for modern
predicted_heights %>% 
  dplyr::select(contains("age2")) %>% 
  rowMeans() -> heights_mod

# Calculate mean percent increase from ancestral to modern
mean((heights_mod - heights_anc)/heights_anc) # 0.003327307
# Calculate 95% credible interval percent increase from ancestral to modern
quantile((heights_mod - heights_anc)/heights_anc, c(0.025, 0.975)) # -0.03802610  0.04696763   


## Collect predicted means for aboveground biomass for CMEM analysis ####
# Collect MCMC samples for all regression coefficients from root-to-shoot model
ggs(agb_cs_out) %>% 
  filter(substr(Parameter, 1, 4) == "beta" | substr(Parameter, 1, 5) == "mu.al") %>% 
  spread(key = Parameter, value = value) -> agb_betas

# Because weight ic_weight and ln_depth were centered, their mean = 0 so they
# don't need to be included in the predicted equation
agb_betas %>% 
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
         frame4.age2.loc2.pred = mu.alpha + `beta[6]` + `beta[7]`) -> predicted_agb

## Collect predicted means for root distribution parameter for CMEM analysis ####
# Collect MCMC samples for all regression coefficients from root-to-shoot model
ggs(beta_cs_out) %>% 
  filter(substr(Parameter, 1, 4) == "beta" | substr(Parameter, 1, 5) == "alpha") %>% 
  spread(key = Parameter, value = value) -> beta_betas

# Because weight ic_weight and ln_depth were centered, their mean = 0 so they
# don't need to be included in the predicted equation
beta_betas %>% 
  mutate(frame1.age1.loc1.pred = alpha + `beta[3]`,
         frame1.age1.loc2.pred = alpha + `beta[3]` + `beta[6]`,
         frame1.age2.loc1.pred = alpha + `beta[3]` + `beta[7]`,
         frame1.age2.loc2.pred = alpha + `beta[3]` + `beta[6]` + `beta[7]`,
         frame2.age1.loc1.pred = alpha + `beta[4]`,
         frame2.age1.loc2.pred = alpha + `beta[4]` + `beta[6]`,
         frame2.age2.loc1.pred = alpha + `beta[4]` + `beta[7]`,
         frame2.age2.loc2.pred = alpha + `beta[4]` + `beta[6]` + `beta[7]`,
         frame3.age1.loc1.pred = alpha + `beta[5]`,
         frame3.age1.loc2.pred = alpha + `beta[5]` + `beta[6]`,
         frame3.age2.loc1.pred = alpha + `beta[5]` + `beta[7]`,
         frame3.age2.loc2.pred = alpha + `beta[5]` + `beta[6]` + `beta[7]`,
         frame4.age1.loc1.pred = alpha,
         frame4.age1.loc2.pred = alpha + `beta[6]`,
         frame4.age2.loc1.pred = alpha + `beta[7]`,
         frame4.age2.loc2.pred = alpha + `beta[6]` + `beta[7]`) -> predicted_beta


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

# Run additive function for each trait (each trait will take ~10 seconds to run)
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
pred_obs_diff <- matrix(NA, nrow = nrow(agb_additive$MonoPredict),
                        ncol = ncol(agb_additive$MonoPredict))

# Loop through polyculture pots (n = 48) to get differences at each iteration.
# This is for figure 2a.
average_difference <- function(trait, additive_samples){
  for (j in 1:48){
    pot_trait_temp <- poly_traits %>%
      filter(pot == j) %>%
      select(trait) %>% as.numeric()
    # This is put in place to make sure that the one pot (39) does not get a
    # beta calculated and it would stand out numerically if it did (i.e. why the
    # temp value is 10000)
    if(trait == "beta" & is.na(pot_trait_temp)){
      pot_trait_temp <- 10000
    }
    # Observed - Predicted (positive values mean that what we observed is greater
    # than the additive expectation, so POSITIVE = FACTILITATION; and negative
    # values mean that what we observed was less than the additive expectation, so
    # NEGATIVE = COMPETITION)
    if(trait == "density"){
      # Predicted values from the poisson regression are on the log-scale
      pred_obs_diff[,j] <- log(pot_trait_temp) - additive_samples$MonoPredict[,j]
    }else{
      pred_obs_diff[,j] <- pot_trait_temp - additive_samples$MonoPredict[,j]
    }
     
  }
  # Then take a mean across all 48 pots for each iteration to get the mean
  # effect (for beta make sure to skip pot 39)
  if(trait == "beta"){
    # Take out pot 39 which has an NA for the observed beta
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
  difference <- matrix(0, nrow = nrow(additive_samples$MonoPredict), ncol = ncol(additive_samples$MonoPredict))
  
  if(trait == "density"){
    # Predicted values from the poisson regression are on the log-scale
    for (i in 1:nrow(additive_samples$MonoPredict)){
      difference[i,] <- log(observed) - additive_samples$MonoPredict[i,]
    }
  }else{
    for (i in 1:nrow(additive_samples$MonoPredict)){
      difference[i,] <- observed - additive_samples$MonoPredict[i,]
    }
  }
  
  # Get average difference across pots for each iteration
  avg_difference <- data.frame(x =colMeans(difference, na.rm = T))
  
  return(avg_difference)
}

# Calculate average difference for each pot for each trait
diffs_agb <- mean_difference_bypot("agb", agb_additive)
diffs_density <- mean_difference_bypot("density", density_additive)
diffs_height<- mean_difference_bypot("mean_tot_height", height_additive)
diffs_width <- mean_difference_bypot("mean_mid_width", width_additive)
diffs_bgb <- mean_difference_bypot("bgb", bgb_additive)
diffs_rs <- mean_difference_bypot("rs", rs_additive)
diffs_beta <- mean_difference_bypot("beta", beta_additive)

# Create a data frame with all average differences and age cohort information
tibble(poly_traits) %>%
  arrange(pot) %>%
  mutate(`aboveground biomass (g)` = diffs_agb$x,
         `stem density` = diffs_density$x,
         `mean stem height (cm)` = diffs_height$x,
         `mean stem width (mm)` = diffs_width$x,
         `belowground biomass (g)` = diffs_bgb$x,
         `root:shoot ratio` = diffs_rs$x,
         `root distribution parameter` = diffs_beta$x) %>%
  # Set the "make-up" of the pot. Pots 1-18 are mixes of four ancestral
  # genotypes, 19-30 are mixes of four descendant genotypes, and 31-48 are mixes
  # of two ancestral and two descendant genotypes.
  mutate(age = case_when(pot %in% 1:18 ~ "ancestral",
                         pot %in% 19:30 ~ "modern",
                         T ~ "mix")) %>%
  dplyr::select(age, `aboveground biomass (g)`, `stem density`, `mean stem height (cm)`, `mean stem width (mm)`,
                `belowground biomass (g)`, `root:shoot ratio`, `root distribution parameter`) -> diffs_by_age

# Save table for later
saveRDS(diffs_by_age, here("outputs/monoculture_polyculture/", "diffs_by_age.rds"))

# Check to see if there are significant differences by age group. The response
# variable here is the scaled difference (difference between obs and pred trait
# values) and the predictor variable is age cohort (ancestral, descendant or
# mixed)
abg_mod <- lm(`aboveground biomass (g)` ~ age, data = diffs_by_age)
anova(abg_mod) 

bgb_mod <- lm(`belowground biomass (g)` ~ age, data = diffs_by_age)
anova(bgb_mod) 

density_mod <- lm(`stem density` ~ age, data = diffs_by_age)
anova(density_mod) 

# Mix is most different here so make it the reference level
diffs_by_age_beta <-
  diffs_by_age %>%
  mutate(age = factor(age)) %>% 
  mutate(age = relevel(age, ref = "mix"))
beta_mod <- lm(`root distribution parameter` ~ age, data = diffs_by_age_beta)
anova(beta_mod)
# Get estimate of mean difference between groups
coef(beta_mod)
# Get confidence intervals around that difference
confint(beta_mod)

rs_mod <- lm(`root:shoot ratio` ~ age, data = diffs_by_age)
anova(rs_mod)
# Get estimate of mean difference between groups
coef(rs_mod)
# Get confidence interval around that difference
confint(rs_mod)

height_mod <- lm(`mean stem height (cm)` ~ age, data = diffs_by_age)
anova(height_mod) 

width_mod <- lm(`mean stem width (mm)` ~ age, data = diffs_by_age)
anova(width_mod) 

## Cohort Marsh Equilibrium Model Simulations - All traits vary ####

## Cohort Marsh Equilibrium Model Simulations

# Parameter estimates for bMax and root:shoot from Blue Genes 2019 experiment
blue_genes <- read_rds(here("supp_data", "blue_genes_subdata.rds"))

# Fit a parabola for the aboveground biomass data from Blue Genes experiment
quad_mod <- lm(agb_scam ~ elevation + I(elevation^2), data = blue_genes)
# Extract quadratic regression coefficients
coefs <- as.numeric(coef(quad_mod))

# Create a function to solve for roots of quadratic formula
quadraticRoots <- function(a, b, c) {
  discriminant <- (b^2) - (4*a*c)
  x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
  x_int_neg <- (-b - sqrt(discriminant)) / (2*a)
  xints <- c(x_int_plus, x_int_neg)
  return(xints)
}

# Extract roots (these are the same as the min and max elevations at which
# biomass can exist)
roots <- quadraticRoots(coefs[3], coefs[2], coefs[1])
zMax_for_sim <- roots[2]
zMin_for_sim <- roots[1]

# Find elevation where peak biomass is given a symmetric parabola
zPeak <- (zMax_for_sim + zMin_for_sim) / 2

# Calculate the predicted biomass at that elevation (bMax)
bMax <- predict(quad_mod, newdata = data.frame(elevation = zPeak))
# Convert to g / cm2 (4 inch diameter pots = 2 inch radius = 5.08 cm radius)
pot_area_cm2 <- pi * 5.08^2
# Get to the scales used in CMEM (g/cm2)
bMax_for_sim <- bMax / pot_area_cm2

# Calculate average root-to-shoot ratio (Note: CMEM assumes static root-to-shoot
# ratio across elevations)

blue_genes %>% 
  # If agb_scam is really low, then root:shoot will be really high just because
  # of the initial propagule's rhizome weight, so we filter out for agb_scam >
  # 0.1
  filter(agb_scam > 0.1) %>% 
  mutate(rs = total_bg / agb_scam) %>% 
  pull(rs) %>% 
  mean() -> rootShoot_for_sim

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
biomass_sum <- get_cv(read_rds(here("outputs/monoculture_models/", "agb_monomodel.rds")))
rs_sum <- get_cv(read_rds(here("outputs/monoculture_models/", "rs_monomodel.rds")))
beta_sum <- get_cv(read_rds(here("outputs/monoculture_models/", "beta_monomodel.rds")))

# Make a vector of the means for each parameter (agb and rs come from blue genes
# and beta comes from our data)
mean_vec <- c(bMax_for_sim, rootShoot_for_sim, beta_sum$mu.alpha)

# Scale biomass to g/cm^2
# Our pots were 6 inches in diameter, 3 inches radius, 7.62 cm radius
pot_area <- 7.62^2*pi
biomass_sum %>% 
  mutate(mu.alpha.scaled = mu.alpha / pot_area,
         sigma.int.scaled = sigma.int / pot_area) -> biomass_sum

##
# All variation - vary due to genotype 
##

# Create covariance matrix if we assume that variance is the same, but means are
# different (conservative assumption)

# Calculate constants that quantify how much larger bMax and rootToShoot are in
# Blue Genes versus this experiment
biomass_scale <- mean_vec[1] / biomass_sum$mu.alpha.scaled
rs_scale <- mean_vec[2] / rs_sum$mu.alpha

# Calculate variances due to genotype given this scaling (following standard
# deviation and variance rules)
var_biomass <- (biomass_sum$sigma.int.scaled * biomass_scale)^2
var_rs <- (rs_sum$sigma.int * rs_scale)^2
# Keep beta the same because the mean is not off of what we might expect
var_beta <- beta_sum$sigma.int^2 

# Calculate covariances between traits based on standard deviations due to
# genotype as well as correlations between traits
covar_biomass_rs <- biomass_sum$sigma.int.scaled * rs_sum$sigma.int * cor(mono_traits$agb, mono_traits$rs)
covar_biomass_beta <- biomass_sum$sigma.int.scaled * beta_sum$sigma.int * cor(mono_traits$agb, mono_traits$beta)
covar_rs_beta <- rs_sum$sigma.int * beta_sum$sigma.int * cor(mono_traits$rs, mono_traits$beta)

# Create covariance matrix
covar_matrix <- matrix(c(var_biomass, covar_biomass_rs, covar_biomass_beta,
                         covar_biomass_rs, var_rs, covar_rs_beta,
                         covar_biomass_beta, covar_rs_beta, var_beta),
                       nrow = 3, byrow = T)

# Take 1000 draws from this multivariate normal distribution
# Set seed to keep this consistent with the manuscript
set.seed(seed)
samples1 <- rmvnorm(1000, mean = mean_vec, covar_matrix)

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
       `maximum rooting depth (cm)` = rooting_depth) -> for_MEM

# Save these samples for later
saveRDS(for_MEM, here("outputs/CMEM_runs/", "traits_for_MEM_simulations.rds"))

# Create storage for MEM model runs
n_runs <- 1000
# Store surface elevation
run_store <- matrix(NA, nrow = n_runs, ncol = n_runs)
# Stores carbon sequestration
carbon_store <- matrix(NA, nrow = n_runs, ncol = n_runs)

# Loop through all of the possible draws (n = 1000) of trait values and push
# through CMEM. Store the elevation between time = 0 (i.e. 2020) and time = 80
# (i.e. 2100) and also calculate the amount of carbon mass at each timestep. (n
# = 1000 iterations takes about ~25 minutes running time)
for (i in 1:n_runs){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0,
                          # Iterate through bMax values
                          bMax = for_MEM$`aboveground biomass (g)`[i], 
                          # Use zMin and zMax (in cm) from Blue Genes experiment
                          zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100,
                          zVegPeak=NA,
                          plantElevationType="orthometric",
                          # Iterate through root-to-shoot ratio
                          rootToShoot = for_MEM$`root:shoot ratio`[i],
                          rootTurnover=0.5,
                          # Iterate through maximum rooting depth
                          rootDepthMax=for_MEM$`maximum rooting depth (cm)`[i],
                          omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  # Calculate amount of carbon at each time step. This follows the %C
  # calculation from Craft et al. 1991
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store[i,]
  print(i)
}

# Calculate average accretion rate for each iteration
init_elev <- 22.6
avg_accretion_rates <- (run_store[,80] - init_elev) / 80

# Calculate average carbon accumulation rate for each iteration
avg_C_accum_rate <- (carbon_store[,80] - carbon_store[,1]) / 80

# Create a data frame of average accretion rates and carbon accumulation rates
# from these simulations for plotting later.
avg_rates <- tibble(avg_acc = avg_accretion_rates,
                    avg_C = avg_C_accum_rate)

write_rds(avg_rates, here("outputs/CMEM_runs", "CMEM_rates_full.rds"))

##
# All variation - vary by cohort mean 
##

# Subset data for just corn and sellman locations
cs_traits <- mono_traits %>% filter(location %in% c('corn', "sellman"))

# Calculate predicted means from cs models for each provenance x cohort
# combination for rs

predicted_rs %>% 
  select(contains("loc1") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_rs_ca
  
predicted_rs %>% 
  select(contains("loc2") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_rs_sa

predicted_rs %>% 
  select(contains("loc1") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_rs_cm

predicted_rs %>% 
  select(contains("loc2") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_rs_sm
  
# Put means in vector: order is corn-ancestral, sellman-ancestral, corn-modern,
# sellman-modern
means_rs <- c(mean_rs_ca, mean_rs_sa, mean_rs_cm, mean_rs_sm)
# Scale based on blue genes means
bg_rs_scale <- rootShoot_for_sim / mean(means_rs)
root_shoot_cohort_forMEM <- means_rs * bg_rs_scale

# Calculate predicted means from cs models for each provenance x cohort
# combination for agb

predicted_agb %>% 
  select(contains("loc1") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_agb_ca

predicted_agb %>% 
  select(contains("loc2") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_agb_sa

predicted_agb %>% 
  select(contains("loc1") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_agb_cm

predicted_agb %>% 
  select(contains("loc2") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_agb_sm

# Put means in vector: order is corn-ancestral, sellman-ancestral, corn-modern,
# sellman-modern and convert to g/cm2
means_agb <- c(mean_agb_ca, mean_agb_sa, mean_agb_cm, mean_agb_sm)/ pot_area

# Scale based on blue genes means
bg_agb_scale <- bMax_for_sim / mean(means_agb)
agb_cohort_forMEM <- means_agb * bg_agb_scale

# Calculate predicted means from cs models for each provenance x cohort
# combination for root distribution parameter (beta)

predicted_beta %>% 
  select(contains("loc1") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_beta_ca

predicted_beta %>% 
  select(contains("loc2") & contains("age1")) %>% 
  rowMeans() %>% mean() -> mean_beta_sa

predicted_beta %>% 
  select(contains("loc1") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_beta_cm

predicted_beta %>% 
  select(contains("loc2") & contains("age2")) %>% 
  rowMeans() %>% mean() -> mean_beta_sm

# Put means in vector: order is corn-ancestral, sellman-ancestral, corn-modern,
# sellman-modern and convert to g/cm2
means_betas <- c(mean_beta_ca, mean_beta_sa, mean_beta_cm, mean_beta_sm)
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

# Run model for each cohort and store surface elevation and carbon accumulation
run_store_cohort <- matrix(NA, nrow = 4, ncol = 100)
carbon_store_cohort <- matrix(NA, nrow = 4, ncol = 100)

for (i in 1:4){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = agb_cohort_forMEM[i], 
                          zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = root_shoot_cohort_forMEM[i],
                          rootTurnover=0.5, rootDepthMax=rooting_depth[i], omDecayRate=0.8,
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

# Calculate CIs around average accretion rate - everything is in the same order
# as before: corn-ancestral, sellman-ancestral, corn-modern, sellman-modern

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

write_rds(cohort_summary, here("outputs/CMEM_runs", "CMEM_rates_full_cohort.rds"))

# Join together data frame of simulations based on genotypic variation and
# cohort simulations
elevation_store_df <- as.data.frame(run_store)
elevation_store_withcohorts <- rbind(elevation_store_df[,1:80], run_store_cohort[,1:80])
colnames(elevation_store_withcohorts) <- 2021:2100

# Reformat data into long format
tibble(elevation_store_withcohorts) %>% 
  mutate(iteration = as.character(1:nrow(elevation_store_withcohorts))) %>% 
  gather(key = year, value = value, `2021`:`2100`) %>% 
  mutate(iteration = case_when(iteration == "1001" ~ "corn-ancestral",
                               iteration == "1002" ~ "sellman-ancestral",
                               iteration == "1003" ~ "corn-modern",
                               iteration == "1004" ~ "sellman-modern",
                               T ~ iteration))  -> CMEM_predictions_belowground

##
# Calculations for in-text
##

# Differences in final marsh elevation due to genotype
CMEM_predictions_belowground %>% 
  filter(year == 2100) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

# Average vertical accretion rates
avg_rates %>% 
  gather(key = state, value = rate, avg_acc:avg_C) %>% 
  group_by(state) %>% 
  summarize(mean = mean(rate),
            lower = quantile(rate, 0.025),
            upper = quantile(rate, 0.975),
            perc_inc = (upper - lower)/lower)

# % diff accretion rate for provenances
cohort_summary %>% 
  group_by(location) %>% 
  summarize(across(where(is.numeric), mean)) -> summary_rates_provenances

# Corn is 1, Sellman is 2
# Vertical accretion rate
(summary_rates_provenances[1, "acc_v"]-summary_rates_provenances[2, "acc_v"])/
  summary_rates_provenances[2, "acc_v"]

# Carbon accumulation rate
(summary_rates_provenances[1, "acc_C"]-summary_rates_provenances[2, "acc_C"])/
  summary_rates_provenances[2, "acc_C"]

# Same for age cohorts
cohort_summary %>% 
  group_by(age) %>% 
  summarize(across(where(is.numeric), mean)) -> summary_rates_cohorts

# Ancestral is 1, Modern is 2
# Vertical accretion rate
(summary_rates_cohorts[1, "acc_v"]-summary_rates_cohorts[2, "acc_v"])/
  summary_rates_cohorts[2, "acc_v"]

# Carbon accumulation rate
(summary_rates_cohorts[1, "acc_C"]-summary_rates_cohorts[2, "acc_C"])/
  summary_rates_cohorts[2, "acc_C"]

## Cohort Marsh Equilibrium Model Simulations - AGB varies only ####

# Take draws from normal distribution with just aboveground biomass varying 
set.seed(seed)
samples_agb <-rnorm(1000, mean = mean_vec[1], sqrt(var_biomass))

# Translate mean beta value to maximum rooting depth
depth_interval <- seq(0, 50, length.out = 1000)
rooting_depth <- NULL

temp_beta <- mean_vec[3]
temp_cumulative <- 1 - temp_beta ^ depth_interval
depth_95 <- which.min(abs(temp_cumulative - 0.95))
rooting_depth <- depth_interval[depth_95]

# Put all data for simulations into a data frame
tibble(`aboveground biomass (g)` = samples_agb,
       `root:shoot ratio` = mean_vec[2],
       `maximum rooting depth (cm)` = rooting_depth) -> for_MEM_agb_only

# Save these samples for later
saveRDS(for_MEM_agb_only, here("outputs/CMEM_runs/", "traits_for_MEM_simulations_agb_only.rds"))

# Create storage for MEM model runs
n_runs <- 1000
run_store_agb_only <- matrix(NA, nrow = n_runs, ncol = 100)
carbon_store_agb_only <- matrix(NA, nrow = n_runs, ncol = 100)

# Run the model for different values of aboveground biomass (n = 1000 iterations
# takes about ~25 minutes running time)
for (i in 1:n_runs){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = for_MEM_agb_only$`aboveground biomass (g)`[i], 
                          zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = for_MEM_agb_only$`root:shoot ratio`[1],
                          rootTurnover=0.5, rootDepthMax=for_MEM_agb_only$`maximum rooting depth (cm)`[1], omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_agb_only[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_agb_only[i,]
  print(i)
}

# Calculate average accretion rates at each iteration
init_elev <- 22.6
avg_accretion_rates_agb_only <- (run_store_agb_only[,80] - init_elev) / 80

# Calculate average carbon accumulation rates at each iteration
avg_C_accum_rate_agb_only <- (carbon_store_agb_only[,80] - carbon_store_agb_only[,1]) / 80

# Allow just agb to vary - vary due to cohort 
run_store_cohort_agb_only <- matrix(NA, nrow = 4, ncol = 100)
carbon_store_cohort_agb_only <- matrix(NA, nrow = 4, ncol = 100)

# For root:shoot and rooting depth, take means across cohorts
for (i in 1:4){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                          initElv=22.6, meanSeaLevel=-1.6,
                          meanHighWaterDatum=12.9, suspendedSediment=3e-05,
                          lunarNodalAmp=0, bMax = agb_cohort_forMEM[i], 
                          zVegMin=zMin_for_sim*100, zVegMax=zMax_for_sim*100, zVegPeak=NA,
                          plantElevationType="orthometric", rootToShoot = mean(root_shoot_cohort_forMEM),
                          rootTurnover=0.5, rootDepthMax=mean(rooting_depth), omDecayRate=0.8,
                          recalcitrantFrac=0.2, captureRate = 2.8)
  run_store_cohort_agb_only[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store_cohort_agb_only[i,]
  print(i)
}

# Calculate average accretion rate
init_elev <- 22.6
avg_accretion_rates_cohort_agb_only <- (run_store_cohort_agb_only[,80] - init_elev) / 80

# Calculate average carbon accumulation rate
avg_C_accum_rate_cohort_agb_only <- (carbon_store_cohort_agb_only[,80] - carbon_store_cohort_agb_only[,1]) / 80

# Bring together information into a data frame
tibble(location = c("corn", "sellman", "corn", "sellman"),
       age = c("ancestral", "ancestral", "modern", "modern"),
       acc_v = avg_accretion_rates_cohort_agb_only * 10,
       acc_C = avg_C_accum_rate_cohort_agb_only * 1e-6 / 1e-8) -> cohort_summary_agb_only

# Store together with the genotype simulations 
elevation_store_df_agb_only <- as.data.frame(run_store_agb_only)
elevation_store_withcohorts_agb_only <- rbind(elevation_store_df_agb_only[,1:80], run_store_cohort_agb_only[,1:80])
colnames(elevation_store_withcohorts_agb_only) <- 2021:2100

# Reformat data into long format
tibble(elevation_store_withcohorts_agb_only) %>% 
  mutate(iteration = as.character(1:nrow(elevation_store_withcohorts_agb_only))) %>% 
  gather(key = year, value = value, `2021`:`2100`) %>% 
  mutate(iteration = case_when(iteration == "1001" ~ "corn-ancestral",
                               iteration == "1002" ~ "sellman-ancestral",
                               iteration == "1003" ~ "corn-modern",
                               iteration == "1004" ~ "sellman-modern",
                               T ~ iteration))  -> CMEM_predictions_agb_only

# Calculate uncertainty at year 2100 for agb only
CMEM_predictions_agb_only %>% 
  filter(year == 2100) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))

# Save output from both simulations for plotting later
write_rds(CMEM_predictions_belowground, here("outputs/CMEM_runs", "CMEM_predictions_full.rds"))
write_rds(CMEM_predictions_agb_only, here("outputs/CMEM_runs", "CMEM_predictions_agb_only.rds"))
