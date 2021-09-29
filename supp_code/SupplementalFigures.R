# Creates all supplemental figures for Belowground Evolution paper

## Preliminaries ####

# Load libraries
library(patchwork);library(raster); library(maps); library(cowplot);
library(ggsn); library(ggmap); library(tidyverse);
library(ggrepel); library(here); library(mvtnorm)

## Read in data
# Derived trait data for all pots
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))
# Read in raw data frame that has seed depths for seed depth to age conversion
for_seed_depths <- read_csv(here("supp_data", "AbovegroundBiomass.csv"))
# Priors needed to predict seed age from seed depth from Vahsen et al.
# germination paper
seed_age_priors <- readRDS(here("supp_data", "SeedAgeCalibrationPriors.rds"))
# Results from monoculture/polyculture analysis
#diffs_by_age <- readRDS("chp1/results/monopoly_diffsbyage.rds")

# Set site colors for mapping and figures below
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Subset data for just monocultures for plots below
mono_traits <- all_traits %>% 
  filter(diversity == "mono")


## Figure S1: map of core locations ####

# Create data frame with GPS points for where cores were taken
locations <- tibble(id = c("Corn Island", "Hog Island", "Kirkpatrick Marsh", "Sellman Creek"),
                    longitude = c(-76.54361, -76.55164, -76.54985, -76.53824),
                    latitude = c(38.87566, 38.87958, 38.87624, 38.89549))

# Create map with scale-bar and N arrow
sbbox <- c(left = -76.555, bottom = 38.87, right = -76.53, top = 38.898)
map <- ggmap(get_stamenmap(sbbox, zoom = 16, maptype = "toner-background")) +
  geom_point(data = locations, aes(x = longitude, y = latitude, color = id), size = 5) +
  scale_color_manual(values = colors) +
  geom_label_repel(data = locations,
                   aes(longitude, latitude,label = id,group = id, color = id),
                   size  = 4, fontface = "bold",
                   box.padding = 0.7, point.padding = 0.5, segment.size = 1.2) +
  scalebar(x.min = -76.555, x.max = -76.532,
           y.min = 38.871,  y.max = 38.872,
           dist = 0.5, transform = TRUE, 
           model = "WGS84", height = 0.3, 
           st.dist = 0.5, dist_unit = "km", st.size = 4, ) +
  xlab("Longitude") +
  ylab("Latitude") + theme_bw() + theme(legend.position = "none")

# Save figure in folder
png(here("figs_tables", "FigS1_map.png"), height = 6, width = 6, res = 300, units = "in")
north2(map, x = 0.28, y = 0.95, scale = 0.06, symbol = 12)
dev.off()

## Figure S2 & Table S1: ages of seeds used ####

# Pull out unique seed depths and order by increasing depth
for_seed_depths %>% 
  dplyr::select(depth_seed) %>% 
  filter(complete.cases(depth_seed)) %>% 
  unique() %>% 
  mutate(depth_seed = as.numeric(depth_seed)) %>% 
  arrange(depth_seed) %>% 
  pull(depth_seed) -> unique_depths

# Get random draws for regression coefficients relating seed depth to seed age
beta_draws <- rmvnorm(1000, seed_age_priors$beta_prior_mean, seed_age_priors$beta_prior_covar)
sigma_draws <- rgamma(1000, seed_age_priors$sigma_prior$alpha, seed_age_priors$sigma_prior$beta)

# Create a matrix to hold predictions of seed age
predicted_age <- matrix(NA, nrow = 1000, ncol = length(unique_depths))

# Loop through to get predicted mean and confidence interval
for(i in 1:1000){
  for (j in 1:length(unique_depths)){
    predicted_age_mean <- beta_draws[i,1] + beta_draws[i,2]*unique_depths[j] + beta_draws[i,3]*unique_depths[j]^2
    predicted_age[i,j] <- rnorm(1, predicted_age_mean, sigma_draws[i])
  }
}

# Reformat from wide to long data
colnames(predicted_age) <- unique_depths
tibble(as.data.frame(predicted_age)) %>% 
  gather(key = depth, value = value, `0`:`16.75`) %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(color_id = ifelse(depth > 5.25, 1, 0))-> pred_age_long

# Get summary statistics for plotting
pred_age_long %>% 
  group_by(depth, color_id) %>% 
  summarize(lower = quantile(value, 0.025),
            median = median(value),
            upper = quantile(value, 0.975)) -> pred_age_summary_stat

# Calculate predicted decade based on seed ages for supplemental table (Table
# S1)
pred_age_summary_stat %>% 
  mutate(pred_year = round(2016 - median, - 1)) -> TableS1

write_csv(TableS1, here("figs_tables", "TableS1.csv"))

# Create figure to show predicted ages by seed depth for the two different age
# cohorts
FigS2 <- pred_age_summary_stat %>% 
  ggplot(aes(x = depth, y = median, group = depth)) +
  geom_pointrange(aes(ymin = lower, ymax = upper, shape = factor(color_id),
                      color = factor(color_id)), size = 1.0, alpha = 0.8) +
  theme_classic() + ylab("predicted seed age (years)") +
  theme(legend.position = "none") +
  coord_flip() +
  scale_x_reverse() +
  scale_color_manual(values = c("gray67", "gray27")) +
  xlab("seed depth (cm)") +
  geom_label(aes(y = 5, x = 2.5), label = "modern", size = 6,
             color = "gray67", label.size = 1) +
  geom_label(aes(y = 75, x = 11.5), label = "ancestral", size = 6,
             color = "gray27", label.size = 1)

png(here("figs_tables", "FigS2_SeedAges.png"), height = 3.5, width = 3.5,
    res = 300, units = "in")
FigS2
dev.off()

## Figure S3: all traits by cohort with raw data ####
tibble(cohort = mono_traits$age,
       provenance = mono_traits$location,
       `aboveground biomass (g)` = mono_traits$agb,
       `stem density` = mono_traits$density,
       `mean stem height (cm)` = mono_traits$mean_tot_height,
       `mean stem width (mm)` = mono_traits$mean_mid_width,
       `belowground biomass (g)` = mono_traits$bgb,
       `root:shoot ratio` = mono_traits$rs,
       `root distribution parameter` = mono_traits$beta) %>% 
  gather(key = trait, value = value, `aboveground biomass (g)`:`root distribution parameter`) %>% 
  mutate(combo = paste(substr(provenance, 1, 1), substr(cohort, 1, 1), sep = "")) %>% 
  mutate(combo = factor(combo, levels = c("ca", "sa", "cm", "sm", "hm", "km"))) %>% 
  mutate(trait = factor(trait, levels = c("aboveground biomass (g)",
                                          "stem density",
                                          "mean stem height (cm)",
                                          "mean stem width (mm)",
                                          "belowground biomass (g)",
                                          "root:shoot ratio",
                                          "root distribution parameter"))) %>% 
  ggplot(aes(x = combo, y = value)) +
  geom_point(aes(color = provenance, shape = cohort),size = 2, alpha = 0.2) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", 
               colour = "black", width = 0.1, 
               position=position_nudge(x = 0.3)) +
  stat_summary(fun = mean, geom = "point", 
               size = 3, aes(color = provenance, shape = cohort),
               position=position_nudge(x = 0.3)) + 
  facet_wrap(~trait, scales = "free_y", nrow = 2) + 
  scale_color_manual(values = colors) +
  ylab("trait value") +
  xlab("") + theme_bw() -> FigS3

png(here("figs_tables", "FigureS3_allTraitsRaw.png"), height = 4, width = 10,
    res = 300, units = "in")
FigS3
dev.off()

## Figure S4: root-to-shoot ratio differs by provenance and cohort ####
mono_traits %>% 
  filter(location %in% c("corn", "sellman")) %>% 
  rename(provenance = location,
         cohort = age) %>% 
  ggplot(aes(x = provenance, y = rs)) +
  geom_boxplot(aes(color = provenance)) +
  geom_jitter(aes(color = provenance, shape = cohort), height = 0, width = 0.1,size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("root-to-shoot ratio") -> rs_location

mono_traits %>% 
  filter(location %in% c("corn", "sellman")) %>% 
  rename(provenance = location,
         cohort = age) %>% 
  ggplot(aes(x = cohort, y = rs)) +
  geom_boxplot() +
  geom_jitter(aes(color = provenance, shape = cohort), height = 0, width = 0.1,size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("") + theme(legend.position = "none") -> rs_age

png(here("figs_tables", "FigS4_rsLocAge.png"), height = 2.7, width = 5.4,
    res = 300, units = "in")
rs_location + rs_age + plot_layout(guides = "collect")
dev.off()

## Figure S5: stem height differs by cohort ####
mono_traits %>% 
  filter(location %in% c("corn", "sellman")) %>% 
  rename(provenance = location,
         cohort = age) %>% 
  ggplot(aes(x = provenance, y = mean_tot_height)) +
  geom_boxplot(aes(color = provenance), outlier.shape = NA) +
  geom_jitter(aes(color = provenance, shape = cohort), height = 0, width = 0.1,
              size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("stem height (cm)") -> height_location

png(here("figs_tables", "FigureS5_heightLoc.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
height_location
dev.off()

## Figure S6: stem width differs by cohort ####
mono_traits %>% 
  filter(location %in% c("corn", "sellman")) %>% 
  rename(provenance = location,
         cohort = age) %>% 
  ggplot(aes(x = cohort, y = mean_mid_width)) +
  geom_boxplot() +
  geom_jitter(aes(color = provenance, shape = cohort), height = 0, width = 0.1,
              size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("stem width (mm)") -> width_age

png(here("figs_tables", "FigureS6_widthAge.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
width_age
dev.off()

# ## Figure S7 & Table SX: monopoly diffs by age ####
# png("chp1/plots/MS_plots/FigureS7_monopoly_byAge.png", height = 4, width = 10, units = "in", res = 300)
# diffs_by_age %>% 
#   gather(key = trait, value = difference, `aboveground biomass (g)`:`root distribution parameter`) %>% 
#   mutate(trait = factor(trait, levels = c("aboveground biomass (g)",
#                                           "stem density",
#                                           "mean stem height (cm)",
#                                           "mean stem width (mm)",
#                                           "belowground biomass (g)",
#                                           "root:shoot ratio",
#                                           "root distribution parameter"))) %>% 
#   mutate(cohort = age) %>% 
#   ggplot(aes(x = cohort, y = difference, pch = cohort)) +
#   geom_boxplot(outlier.shape = NA)+
#   geom_jitter(height = 0, width = 0.2, alpha = 0.2, size = 3) +
#   facet_wrap(~trait, scales = "free_y", nrow = 2) +
#   scale_shape_manual(values = c(16,8,17)) +
#   ylab("scaled difference")
# dev.off()
# # Also test for differences across age cohorts
# anova(lm(`aboveground biomass (g)` ~ age, data = diffs_by_age)) # ns
# anova(lm(`stem density` ~ age, data = diffs_by_age)) # ns
# anova(lm(`mean stem height (cm)` ~ age, data = diffs_by_age)) # ns
# anova(lm(`mean stem width (mm)` ~ age, data = diffs_by_age)) # ns
# anova(lm(`belowground biomass (g)` ~ age, data = diffs_by_age)) # .
# anova(lm(`root:shoot ratio` ~ age, data = diffs_by_age)) # *
# anova(lm(`root distribution parameter` ~ age, data = diffs_by_age)) # ns
# 
# 
# 
# ## Figure S8: trait space for MEM simulations ####
# for_MEM1 <- readRDS("chp1/results/randomdrawsMEM.rds")
# 
# png("chp1/plots/MS_plots/FigureS8_randomdrawsMEM.png", width = 6.7, height = 5.7, res = 300, units = "in")
# GGally::ggpairs(for_MEM1, lower = list(continuous = wrap("smooth", alpha = 0.3, size=1)))
# dev.off()
