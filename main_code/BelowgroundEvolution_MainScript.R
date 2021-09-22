# Vahsen et al. Belowground Evolution main script

# This script runs the main analyses and sources additional functions within the
# main_code folder

## Load libraries ####
library(tidyverse); library(ggmcmc); library(rjags)

## Read in data ####

# Pot-level trait data for (generalized) linear mixed models
pot_traits <- read_csv("data/BelowgroundEvolution_PotLevel.csv")

# Source JAGS functions to fit (generalized) linear mixed models
