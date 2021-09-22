# Vahsen et al. Belowground Evolution main script

# This script runs the main analyses and sources additional functions within the
# main_code folder

## Load libraries ####
library(tidyverse); library(ggmcmc); library(rjags); library(lme4)

## Read in data ####

# Pot-level trait data for (generalized) linear mixed models
pot_traits <- read_csv("data/BelowgroundEvolution_PotLevel.csv")

## Data formatting ####
# Set appropriate factors
pot_traits %>% 
  mutate(pot = factor(pot),
         frame = factor(frame),
         diversity = factor(diversity),
         age = factor(age),
         provenance = factor(provenance),
         genotype = factor(genotype)) -> pot_traits
# Subset data into all monoculture data and all polyculture data
pot_traits %>% 
  filter(diversity == "mono") -> mono_traits
pot_traits %>% 
  filter(diversity == "poly") -> poly_traits


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
saveRDS(agb_out, "outputs/agb_monomodel.rds")
saveRDS(rs_out, "outputs/rs_monomodel.rds")
saveRDS(beta_out, "outputs/beta_monomodel.rds")
saveRDS(width_out, "outputs/width_monomodel.rds")
saveRDS(height_out, "outputs/height_monomodel.rds")
saveRDS(rs_out, "outputs/rs_monomodel.rds")
saveRDS(bgb_out, "outputs/bgb_monomodel.rds")

## Fit all trait models for age + provenance additive effects ####

# Subset data to just Corn Island and Sellman Creek because those provenances
# have multiple reps within each age cohort (only drops 6 reps)
mono_traits %>% 
  filter(provenance %in% c("corn", "sellman")) -> cs_traits

# Previous models fits showed non-significant age by provenance interactions for
# all trait models. So here we fit all models with additive effects of age and
# provenance 

# Create new model template
model_template_cs <- lmer(agb ~ ln_depth + ic_weight + frame +
                            provenance + age + (1|genotype), data = cs_traits)

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

