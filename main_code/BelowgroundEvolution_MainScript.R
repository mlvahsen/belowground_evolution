# Vahsen et al. Belowground Evolution main script

# This script runs the main analyses and sources additional functions within the
# main_code folder

## Load libraries ####
library(tidyverse); library(ggmcmc); library(rjags); library(lme4)

## Read in data ####

# Pot-level trait data for (generalized) linear mixed models
pot_traits <- read_csv("data/BelowgroundEvolution_PotLevel.csv")
# Initial condition data of which propagule were in each pot
ic_stems <- read_csv("data/InitialPropagules.csv")

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

# Save as R data file objects to use later for plotting and subsequent analyses
saveRDS(agb_cs_out, "outputs/agb_csmodel.rds")
saveRDS(rs_cs_out, "outputs/rs_csmodel.rds")
saveRDS(beta_cs_out, "outputs/beta_csmodel.rds")
saveRDS(width_cs_out, "outputs/width_csmodel.rds")
saveRDS(height_cs_out, "outputs/height_csmodel.rds")
saveRDS(rs_cs_out, "outputs/rs_csmodel.rds")
saveRDS(bgb_cs_out, "outputs/bgb_csmodel.rds")
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
mean(rs_corn - rs_sellman) # 0.1023549
# Calculate 95% credible interval of difference in means
quantile(rs_corn - rs_sellman, c(0.025, 0.975)) # 0.01852034 0.18912246
# Calculate mean percent increase from Sellman to Corn
mean(rs_corn / rs_sellman - 1) # 0.1458952
# Calculate 95% credible interval percent increase from Sellman to Corn
quantile(rs_corn/rs_sellman - 1, c(0.025, 0.975)) # 0.02368095 0.28810277 

# Calculate average predicted root-to-shoot for ancestral cohort
predicted_rs %>% 
  select(contains("age1")) %>% 
  rowMeans() -> rs_ancestral

# Calculate average predicted root-to-shoot for modern cohort
predicted_rs %>% 
  select(contains("age2")) %>% 
  rowMeans() -> rs_modern

# Calculate mean percent decrease from ancestral to modern
mean((rs_ancestral - rs_modern)/rs_ancestral) # 0.07282689
# Calculate 95% quantile for percent decrease from ancestral to modern
quantile((rs_ancestral - rs_modern)/rs_ancestral, c(0.025, 0.975)) # -0.02792931  0.17022620 

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
mean((widths_ancestral - widths_modern)/widths_ancestral) # 0.05637774
# Calculate 95% credible interval of percent decrease from ancestral to modern
quantile((widths_ancestral - widths_modern)/widths_ancestral, c(0.025, 0.975)) # -0.03070937  0.14518776   

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
mean((heights_sellman - heights_corn)/heights_corn) # 0.02920073
# Calculate 95% credible interval percent increase from Corn to Sellman
quantile((heights_sellman - heights_corn)/heights_corn, c(0.025, 0.975)) # -0.01567764  0.07571127   

## Monoculture vs polyculture analysis ####
## Monoculture vs Polyculture ####
additive_predict <- function(monoculture_ggs){
  # Set up information about polyculture pots (initial weights of propagules,
  # frame and ln_depth). Not actually fitting a model here, just getting
  # information from the model matrix
  mod_poly <- lm(agb ~ ic_weight + ln_depth + frame, data = poly_traits) 
  
  # Get model matrices from linear models (remove intercept)
  M_poly <- as.matrix(model.matrix(mod_poly)[ , -1])
  
  # Scale continuous covariates (using means and sd from in monoculture models)
  M_poly_sc <- M_poly
  ic_weight_Mono_mean <- mean(mono_traits$ic_weight)
  ic_weight_Mono_sd <- sd(mono_traits$ic_weight)
  ln_depth_Mono_mean <- mean(mono_traits$ln_depth)
  ln_depth_Mono_sd <- sd(mono_traits$ln_depth)
  M_poly_sc[,"ic_weight"] <- (M_poly[,"ic_weight"] - ic_weight_Mono_mean)/ic_weight_Mono_sd
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
    
    ic_stems %>% 
      filter(pot == i) %>% 
      mutate(weight_4 = Weight * 4) -> pot_mono
    
    # Pull covariate information
    M_poly_sc[i,2:5] -> pot_poly
    
    # Make data easy to grab from jags_mono object
    monoculture_ggs %>%
      filter(substr(Parameter, 1, 5) == "alpha") %>% 
      spread(key = Parameter, value = value) %>% 
      dplyr::select(contains("alpha")) -> alphas_ggs
    
    alphas_df <- as.data.frame(alphas_ggs)
    colnames(alphas_df) <- levels(mono_traits$genotype)
    
    # Get regression coefficient for ln_depth
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
      filter(Parameter == "beta[2]") %>% 
      pull(value) -> beta_frame2
    # Get last frame regression coefficient
    monoculture_ggs %>% 
      filter(Parameter == "beta[3]") %>% 
      pull(value) -> beta_frame3
    
    
    # Create regression with frame, ic_weight, and depth info from the
    # polyculture and use regression coefficients from the monoculture model
    for (j in 1:4){
      genotype_name <- pull(pot_mono[j, "Id"])
      predict[,j] <- alphas_df[, genotype_name] + 
        beta_icweight * pull(pot_mono[j, "Weight"]) + beta_lndepth * pot_poly[1] +
        beta_frame1 * pot_poly[2] + beta_frame2 * pot_poly[3] + beta_frame3 * pot_poly[4]
    }
    
    # Then take the row means to get the predicted biomass if this was all additive
    MonoPredict[,i] <- rowMeans(predict, na.rm = T)
  }
  return(list(MonoPredict = MonoPredict))
}

biomass_additive <- additive_predict(ggs(agb_out))
bgb_additive <- additive_predict(ggs(bgb_out))
height_additive <- additive_predict(ggs(height_out))
width_additive <- additive_predict(ggs(width_out))
rs_additive <- additive_predict(ggs(rs_out))
density_additive <- additive_predict(ggs(density_out))
beta_additive <- additive_predict(ggs(beta_out))

saveRDS(biomass_additive, "chp1/results/monopoly_biomass_add.rds")
saveRDS(bgb_additive, "chp1/results/monopoly_bgb_add.rds")
saveRDS(height_additive, "chp1/results/monopoly_height_add.rds")
saveRDS(width_additive, "chp1/results/monopoly_width_add.rds")
saveRDS(rs_additive, "chp1/results/monopoly_rs_add.rds")
saveRDS(density_additive, "chp1/results/monopoly_density_add.rds")
saveRDS(beta_additive, "chp1/results/monopoly_beta_add.rds")

mean_difference_bypot <- function(trait, additive_samples){
  # Predicted values for each polyculture pot (row = iterations, col = pots)
  observed <- poly_all[,trait]
  
  # Create matrix to hold differences between observed and predicted for each pot
  # at each iteration
  difference <- matrix(0, nrow = 4000, ncol = 48)
  
  if(trait == "density"){
    for (i in 1:4000){
      difference[i,] <- log(observed) - additive_samples$MonoPredict[i,]
    }
  }else{
    for (i in 1:4000){
      difference[i,] <- observed - additive_samples$MonoPredict[i,]
    }
  }
  
  # Get average difference across pots for each iteration
  avg_difference <- data.frame(x =colMeans(difference/mean(observed, na.rm = T), na.rm = T))
  
  return(avg_difference)
}

# Check for effects of age on difference for each trait
diffs_biomass <- mean_difference_bypot("total_biomass", biomass_additive)
diffs_density <- mean_difference_bypot("density", density_additive)
diffs_height<- mean_difference_bypot("mean_tot_height", height_additive)
diffs_width <- mean_difference_bypot("mean_mid_width", width_additive)
diffs_bgb <- mean_difference_bypot("total_bgb", bgb_additive)
diffs_rs <- mean_difference_bypot("rs", rs_additive)
diffs_beta <- mean_difference_bypot("beta", beta_additive)


tibble(poly_all) %>%
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

saveRDS(diffs_by_age, "chp1/results/monopoly_diffsbyage.rds")

# WHOAAA -- looks like bg biomass is going in the direction we would predict!!
# Competition is not as strong in modern vs ancestral 

abg_mod <- lm(`ag biomass` ~ age, data = diffs_by_age)
anova(abg_mod)

bgb_mod <- lm(`bg biomass` ~ age, data = diffs_by_age)
anova(bgb_mod) # .

density_mod <- lm(density ~ age, data = diffs_by_age)
anova(density_mod)

beta_mod <- lm(`root parameter` ~ age, data = diffs_by_age)
anova(beta_mod)

rs_mod <- lm(`root:shoot` ~ age, data = diffs_by_age)
anova(rs_mod) # *

height_mod <- lm(`stem height` ~ age, data = diffs_by_age)
anova(height_mod)

width_mod <- lm(`stem width` ~ age, data = diffs_by_age)
anova(width_mod)

# Pull out bgb and root:shoot to show systematic differences


