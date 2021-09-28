# Bring in all trait data to prepare for (generalized) linear models,
# monoculture vs polyculture analysis, and all plots

## Preliminaries ####

## Load libraries
library(tidyverse); library(here); library(lme4)

# Make %notin% operator
`%notin%` <- Negate(`%in%`)

# Set iterations for fitting Bayesian models
n.adapt <- 50000
n.iter <- n.adapt*4
thin <- n.adapt/500

## Read in data
# Source code that makes adjustment for stem weights for final belowground
# biomass
source(here("supp_code", "AssembleBelowgroundData.R")) # Warning is ok (MLV checked)
# Read in initial propagule weights
source(here("supp_code", "GCREW_InitialConditions.R"))
# Read in aboveground biomass data
agb <- read_csv(here("supp_data", "AbovegroundBiomass.csv"))
# Read in census data
census <- read_csv(here("supp_data", "FinalCensus.csv"))

## Adjust agb based on separating pieces of stems from sieved belowground
## biomass
agb %>% 
  arrange(pot) %>% 
  mutate(agb = total_biomass + arrange(agb_adjust, pot)) -> agb

# Missing weight for allometric stem 2 for pot 11 so use allometric data to
# calculate what the weight should be.

# Allometric calculation
census %>% 
  filter(complete.cases(allomet) & allomet != "A") %>% 
  mutate(allomet_stemID = paste(pot, allomet, sep = "_")) %>% 
  select(allomet_stemID, tot_height, mid_width) %>% 
  filter(allomet_stemID != "11_2") -> allomets

agb %>% 
  select(pot, allom_1:allom_5) %>% 
  gather(key = stem, value = weight, allom_1:allom_5) %>% 
  mutate(allomet_stemID = paste(pot, parse_number(stem), sep = "_")) %>% 
  filter(allomet_stemID != "11_2") -> allomets_weight

merge(allomets, allomets_weight) -> for_allomet_calc

for_allomet_calc %>% 
  # square-root transform is pretty good
  ggplot(aes(x = tot_height, y = sqrt(weight))) +
  geom_point() 

# Fit model to predict biomass for pot 11 allomet 2
allomet_model <- lmer(sqrt(weight) ~ tot_height * mid_width + (1|pot), data = for_allomet_calc)
plot(allomet_model)
qqnorm(resid(allomet_model))

new_data <- census %>%
  filter(allomet == 2 & pot == 11) %>%
  select(pot, tot_height, mid_width) %>% 
  as.data.frame()

to_add_pot11_agb <- predict(allomet_model, new_data)^2

# Adjust agb for pot 11 by adding prediction based on linear model
agb %>% 
  mutate(total_biomass = case_when(pot == 11 ~ total_biomass + to_add_pot11_agb,
                                   T ~ total_biomass)) -> agb

# Get total initial weight for each pot by adding weights of four propagules
ic_data %>% 
  group_by(pot) %>% 
  dplyr::summarize(ic_weight = sum(weight)) -> initial

# Link up initial conditions data with abg data
agb <- merge(agb, initial, by = "pot") %>% 
  mutate(density = density_orig_live) %>% 
  dplyr::select(pot, frame, old_genotype, location, age,
                ic_weight, depth, agb = total_biomass, density) %>% 
  mutate(genotype = old_genotype) %>% 
  mutate(ln_depth = log(depth))

# Deal with census data
# Remove dead stems and tidy
census %>% 
  filter(dead == 0) %>% 
  dplyr::select(pot, tot_height, mid_width, diversity, bas_width) %>%  
  group_by(pot) %>% 
  dplyr::summarize(mean_mid_width = mean(mid_width),
                   mean_bas_width = mean(bas_width),
                   mean_tot_height = mean(tot_height),
                   max_tot_height = max(tot_height)) -> live_census_sum


# Merge together mono and live_census_sum
census_all <- merge(agb, live_census_sum, by = "pot") 

# Rename bgb for downstream formatting
bgb <- bgb_adjusted_full

# Calculate total belowground biomass and root-to-shoot ratio and merge with
# rest of data
bgb %>% 
  group_by(pot) %>% 
  summarize(bgb = sum(biomass_nostems)) %>% 
  arrange(pot) %>% 
  pull(bgb) -> bgb_vector

census_all %>% 
  arrange(pot) %>% 
  mutate(bgb = bgb_vector,
         rs = bgb / agb) -> census_all

## Start root depth distribution calculations
# Make summary of total biomass and number of layers for each pot
bgb %>%
  # Remove pot 39 because it had half power-washed and half not
  filter(pot != 39) %>% 
  group_by(pot) %>%  
  summarize(total = sum(biomass),
            n = length(biomass)) -> totals

# Make sure data is sorted by pot and rooting depth
bgb %>%
  filter(pot != 39) %>% 
  arrange(pot, depth_roots) %>% 
  # Create repeated entries of the total belowground biomass in order to
  # calculate proportions
  mutate(total = rep(totals$total, totals$n)) %>% 
  # Calculate proportions
  mutate(prop = biomass / total) %>% 
  # Group by pot so that we can calculate cumulative proportion across depth
  group_by(pot) %>% 
  mutate(cum_sum = cumsum(prop)) %>% 
  ungroup() -> bgb_beta_analysis

# Create an additional set of data such that at depth 0, there is no biomass
bgb_beta_analysis %>% 
  dplyr::select(frame, position, pot, age, diversity, genotype, old_genotype, location, depth_seed) %>% 
  unique() -> for_adding_zeros 

# Here depth is -5, but will get changed to 0 in the next step
for_adding_zeros$depth_roots <- -5
for_adding_zeros$biomass <- 0
for_adding_zeros$cum_sum <- 0

# Order columns the same way as the for_adding_zeros dataframe
bgb_beta_analysis %>% 
  dplyr::select(frame, position, pot, age,
                diversity, genotype, old_genotype,
                location, depth_seed, depth_roots, biomass, cum_sum)->bgb_beta_analysis

bgb_beta_analysis_zeros <- rbind(for_adding_zeros, bgb_beta_analysis) %>% 
  mutate(depth_roots_adjust = depth_roots + 5) %>% 
  arrange(pot, depth_roots)

## Run non-linear regression on cumulative probability for each pot

# Create output vector to hold information
beta <- NULL

for (i in 1:length(unique(for_adding_zeros$pot))){
  sub_data <- bgb_beta_analysis_zeros %>% filter(pot == unique(for_adding_zeros$pot)[i])
  mod_temp <- nls(cum_sum ~ 1 - beta ^ (depth_roots_adjust), data = sub_data, start = list(beta = 1))
  beta[i] <- coef(mod_temp)[1]
}

for_adding_zeros$beta <- beta

# Merge together with the rest of the data
for_adding_zeros %>% 
  select(pot, beta) %>% 
  add_row(pot = 39, beta = NA) %>% 
  arrange(pot) -> betas_to_merge

merge(census_all, betas_to_merge, by = "pot") -> census_all

# Create column for diversity 
census_all %>% 
  mutate(diversity = case_when(pot %in% 1:48 ~ "poly",
         T ~ "mono")) %>% 
  relocate(diversity, .after = frame) -> trait_data_derived

# Now we have full dataset for monocultures and polycultures
write_csv(trait_data_derived, here("data", "CompiledTraitData.csv"))
