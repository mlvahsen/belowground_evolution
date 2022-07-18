# Creates all supplemental figures for Belowground Evolution paper

## Preliminaries ####

# Load libraries
library(patchwork);library(raster); library(maps); library(cowplot);
library(ggsn); library(ggmap); library(tidyverse);
library(ggrepel); library(here); library(mvtnorm); library(GGally);
library(ggmcmc); library(ggExtra)

## Read in data
# Derived trait data for all pots
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))
# Read in raw data frame that has seed depths for seed depth to age conversion
for_seed_depths <- read_csv(here("supp_data", "AbovegroundBiomass.csv"))
# Priors needed to predict seed age from seed depth from Vahsen et al.
# germination paper
seed_age_priors <- readRDS(here("supp_data", "SeedAgeCalibrationPriors.rds"))
# Results from monoculture/polyculture analysis
diffs_by_age <- readRDS(here("outputs/monoculture_polyculture/", "diffs_by_age.rds"))
# Blue genes data to inform mean adjustments for CMEM simulation
blue_genes <- read_rds("supp_data/blue_genes_subdata.rds")

# Set site colors for mapping and figures below
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Subset data for just monocultures for plots below
mono_traits <- all_traits %>% 
  filter(diversity == "mono")

# Read in all monoculture trait models
agb_out <- read_rds(here("outputs/monoculture_models", "agb_monomodel.rds"))
bgb_out <- read_rds(here("outputs/monoculture_models", "bgb_monomodel.rds"))
rs_out <- read_rds(here("outputs/monoculture_models", "rs_monomodel.rds"))
height_out <- read_rds(here("outputs/monoculture_models", "height_monomodel.rds"))
width_out <- read_rds(here("outputs/monoculture_models", "width_monomodel.rds"))
density_out <- read_rds(here("outputs/monoculture_models", "density_monomodel.rds"))
beta_out <- read_rds(here("outputs/monoculture_models", "beta_monomodel.rds"))

## Figure S1: map of core locations ####

# Create data frame with GPS points for where cores were taken
locations <- tibble(id = c("Corn Island", "Hog Island", "Kirkpatrick Marsh", "Sellman Creek"),
                    longitude = c(-76.54361, -76.55164, -76.54985, -76.53824),
                    latitude = c(38.87566, 38.87958, 38.87624, 38.89549))

# Create map with scale-bar and N arrow
sbbox <- c(left = -76.555, bottom = 38.87, right = -76.53, top = 38.898)
map <- ggmap(get_stamenmap(sbbox, zoom = 16, maptype = "terrain-background")) +
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
           st.dist = 0.5, dist_unit = "km", st.size = 4) +
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
set.seed(1234)
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
  mutate(pred_year = round(2016 - median)) -> TableS1

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
  geom_label(aes(y = 10, x = 2.5), label = "descendant", size = 6,
             color = "gray67", label.size = 1) +
  geom_label(aes(y = 75, x = 11.5), label = "ancestral", size = 6,
             color = "gray27", label.size = 1)

png(here("figs_tables", "FigS2_SeedAges.png"), height = 3.5, width = 3.5,
    res = 300, units = "in")
FigS2
dev.off()

## Figure S3: random intercepts of monoculture models ####
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

mono_traits %>% 
  rename(`age cohort` = age,
         provenance = location) -> mono_traits

# Calculate ICCs
get_icc <- function(coda_object){
  ggs(coda_object, family = "sigma") %>% 
    spread(key = Parameter, value = value) %>% 
    mutate(icc = sigma.int^2 / (sigma.int^2 + sigma.res^2)) %>% 
    summarize(mean = sprintf('%.1f', mean(icc)*100),
              lower.95 = sprintf('%.1f', quantile(icc, 0.025)*100),
              upper.95 = sprintf('%.1f', quantile(icc, 0.975)*100)) -> icc_summary
  return(icc_summary)
}
get_icc_pois <- function(coda_object){
  ggs(coda_object) %>% 
    spread(key = Parameter, value = value) %>%
    mutate(sigma2_d = log(1 + 1/mean(mono_traits$density))) %>% 
    mutate(icc = sigma.int^2 / (sigma.int^2 + sigma2_d)) %>% 
    summarize(mean = sprintf('%.1f', mean(icc)*100),
              lower.95 = sprintf('%.1f', quantile(icc, 0.025)*100),
              upper.95 = sprintf('%.1f', quantile(icc, 0.975)*100)) -> icc_summary_density
  return(icc_summary_density)
}

# Plotting function
make_fig1_panel <- function(data, coda_object, trait, label_y, legend, xlab){
  data %>% 
    mutate(`age cohort` = case_when(`age cohort` == "modern" ~ "descendant",
                                    `age cohort` == "ancestral" ~ "ancestral")) %>% 
    mutate(genotype_name = genotype) %>% 
    mutate(genotype_code = as.numeric(as.factor(genotype))) %>% 
    group_by(genotype_code, genotype_name, provenance, `age cohort`) %>% 
    dplyr::summarize(n = length(`age cohort`)) %>% 
    mutate(Parameter = paste("alpha[", genotype_code, "]", sep = "")) %>% 
    mutate(Coefficient = "Intercept") %>% 
    mutate(Label = genotype_name) -> labels
  
  merge(ggs(coda_object) %>% filter(substr(Parameter, 1, 5) == "alpha"),
        labels, by.x = "Parameter") -> out
  
  out %>% 
    mutate(new_Label = case_when(Label == "C1BL" ~ "ca1",
                                 Label == "C1BN" ~ "ca2",
                                 Label == "C2BL" ~ "ca3",
                                 Label == "C2BM" ~ "ca4",
                                 Label == "C4BM" ~ "ca5",
                                 Label == "S1CDI" ~ "sa1",
                                 Label == "S1CSJ" ~ "sa2",
                                 Label == "S1CSK" ~ "sa3",
                                 Label == "C1BP" ~ "cd1",
                                 Label == "C1BR" ~ "cd2",
                                 Label == "C2BT" ~ "cd3",
                                 Label == "C3AS" ~ "cd4",
                                 Label == "S1ADN" ~ "sd1",
                                 Label == "S1CDP" ~ "sd2",
                                 Label == "H1A1P" ~ "hd1",
                                 Label == "KM1B2P" ~ "kd1")) -> out
  
  out$new_Label <- factor(out$new_Label , levels=c("ca1", "ca2", "ca3", "ca4","ca5",
                                                   "sa1", "sa2", "sa3",
                                                   "cd1", "cd2", "cd3", "cd4",
                                                   "sd1", "sd2", "hd1", "kd1"))
  
  if(trait == "stem density"){
    out %>% 
      mutate(value = exp(value)) -> out
  }
  
  # Get mean value to plot as hline
  ggs(coda_object) %>% 
    filter(Parameter == "mu.alpha") %>% 
    summarise(mean = mean(value)) %>% pull(mean) -> mean_value
  
  if(trait == "stem density"){
    mean_value <- exp(mean_value)
  }
  
  if(trait == "stem density"){
    icc <- get_icc_pois(coda_object)
  }else{
    icc <- get_icc(coda_object)
  }
  
  
  
  ggplot(out) +
    geom_pointrange(mapping = aes(x = new_Label, y = value, color = provenance, shape = `age cohort`),
                    stat = "summary",
                    fun.min = function(z) {quantile(z,0.025)},
                    fun.max = function(z) {quantile(z,0.975)},
                    fun = mean,
                    size = 0.8) +
    geom_vline(aes(xintercept = 8.5)) +
    geom_hline(aes(yintercept = mean_value), color = "gray47", linetype = "dashed") +
    scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
    xlab("genotype") + ylab(trait) +
    theme_bw() +
    theme(legend.position = legend, legend.box = "vertical",
          legend.title=element_text(size=14, face = "bold"), 
          legend.text=element_text(size=12)) + 
    xlab(xlab) +
    annotate("text", x = 13.5, y = label_y, label = paste("ICC = ", icc$mean, " (",
                                                             icc$lower.95, ",", icc$upper.95, ")", sep = "")) 
    
  
  
  
  
}

# Make plots
fig1_bgb <- make_fig1_panel(mono_traits, bgb_out, "belowground biomass (g)", 8.5, "top", "")
fig1_agb <- make_fig1_panel(mono_traits, agb_out, "aboveground biomass (g)", 9, "none", "")
fig1_density <- make_fig1_panel(mono_traits, density_out, "stem density", 55, "none", "")
fig1_width <- make_fig1_panel(mono_traits, width_out, "mean stem width (mm)", 3.4, "none", "")
fig1_height <- make_fig1_panel(mono_traits, height_out, "mean stem height (cm)", 47.5, "none", "")
fig1_rs <- make_fig1_panel(mono_traits, rs_out, "root:shoot ratio", 1.0, "none", "genotype")
fig1_beta <- make_fig1_panel(mono_traits, beta_out, trait = "root distribution parameter", 0.93, "none", "genotype")

# Pull out legend
legend <- get_legend(fig1_bgb)
fig1_bgb_nolegend <- fig1_bgb + theme(legend.position = "none")

# Put together LHS
left_panel <- cowplot::plot_grid(fig1_agb, 
                                 fig1_density, 
                                 fig1_height,
                                 fig1_width, nrow = 4,
                                 labels = c("a", "c", "e", "g"),
                                 label_x = 0.12, align = "v")

# Put together RHS except legend
right_panela <- cowplot::plot_grid(fig1_bgb_nolegend, 
                                   fig1_rs, 
                                   fig1_beta,
                                   nrow = 3,
                                   labels = c("b", "d", "f"),
                                   label_x = 0.12, align = "v")

# Add legend to RHS
right_panel<- cowplot::plot_grid(right_panela, legend, nrow = 2, rel_heights = c(3,1))

png(here("figs_tables", "FigS3.png"), height = 8.6, width = 10, res = 300, units = "in")
left_panel + right_panel
dev.off() 

## Figure S4: monopoly diffs by age ####
diffs_by_age %>%
  gather(key = trait, value = difference, `aboveground biomass (g)`:`root distribution parameter`) %>%
  mutate(trait = factor(trait, levels = c("aboveground biomass (g)",
                                          "stem density",
                                          "mean stem height (cm)",
                                          "mean stem width (mm)",
                                          "belowground biomass (g)",
                                          "root:shoot ratio",
                                          "root distribution parameter"))) %>%
  mutate(cohort = age) %>%
  mutate(cohort = case_when(cohort == "ancestral" ~ "ancestral",
                            cohort == "modern" ~ "descendant",
                            T ~ "mix")) %>% 
  mutate(cohort = factor(cohort, levels = c("ancestral", "mix", "descendant"))) %>% 
  ggplot(aes(x = cohort, y = difference, pch = cohort)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.2, size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  facet_wrap(~trait, scales = "free_y", nrow = 2) +
  scale_shape_manual(values = c(16,8,17)) +
  ylab("scaled difference") + theme_classic() +
  theme(legend.position = c(1, 0.1),
        legend.justification = c(1.6, 0)) -> fig_S4

png(here("figs_tables","FigureS4_monopoly_byAge.png"), height = 4, width = 10, units = "in", res = 300)
fig_S4
dev.off()

# Also test for differences across age cohorts
anova(lm(`aboveground biomass (g)` ~ age, data = diffs_by_age)) # ns
anova(lm(`stem density` ~ age, data = diffs_by_age)) # ns
anova(lm(`mean stem height (cm)` ~ age, data = diffs_by_age)) # ns
anova(lm(`mean stem width (mm)` ~ age, data = diffs_by_age)) # ns
anova(lm(`belowground biomass (g)` ~ age, data = diffs_by_age)) # ns
anova(lm(`root:shoot ratio` ~ age, data = diffs_by_age)) # ns
anova(lm(`root distribution parameter` ~ age, data = diffs_by_age)) # .

## Figure S5: root-to-shoot ratio differs by provenance and cohort ####
mono_traits %>% 
  filter(provenance %in% c("corn", "sellman")) %>% 
  mutate(`age cohort` = case_when(`age cohort` == "ancestral" ~ "ancestral",
                                  T ~ "descendant")) %>% 
  ggplot(aes(x = provenance, y = rs)) +
  geom_boxplot(aes(color = provenance)) +
  geom_jitter(aes(color = provenance, shape = `age cohort`), height = 0, width = 0.1,size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("root-to-shoot ratio") -> rs_location

mono_traits %>% 
  filter(provenance %in% c("corn", "sellman")) %>% 
  mutate(`age cohort` = case_when(`age cohort` == "ancestral" ~ "ancestral",
                                  T ~ "descendant")) %>%
  ggplot(aes(x = `age cohort`, y = rs)) +
  geom_boxplot() +
  geom_jitter(aes(color = provenance, shape = `age cohort`), height = 0, width = 0.1,size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("") + theme(legend.position = "none") -> rs_age

png(here("figs_tables", "FigS5_rsLocAge.png"), height = 2.7, width = 5.4,
    res = 300, units = "in")
rs_location + rs_age + plot_layout(guides = "collect")
dev.off()

## Figure S6: stem width differs by cohort ####
mono_traits %>% 
  filter(provenance %in% c("corn", "sellman")) %>% 
  mutate(`age cohort` = case_when(`age cohort` == "ancestral" ~ "ancestral",
                                  T ~ "descendant")) %>%
  ggplot(aes(x = `age cohort`, y = mean_mid_width)) +
  geom_boxplot() +
  geom_jitter(aes(color = provenance, shape = `age cohort`), height = 0, width = 0.1,
              size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("stem width (mm)") -> width_age

png(here("figs_tables", "FigureS6_widthAge.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
width_age
dev.off()

## Figure S7: all traits by cohort with raw data ####
tibble(cohort = mono_traits$`age cohort`,
       provenance = mono_traits$provenance,
       `aboveground biomass (g)` = mono_traits$agb,
       `stem density` = mono_traits$density,
       `mean stem height (cm)` = mono_traits$mean_tot_height,
       `mean stem width (mm)` = mono_traits$mean_mid_width,
       `belowground biomass (g)` = mono_traits$bgb,
       `root:shoot ratio` = mono_traits$rs,
       `root distribution parameter` = mono_traits$beta) %>% 
  gather(key = trait, value = value, `aboveground biomass (g)`:`root distribution parameter`) %>% 
  mutate(cohort = case_when(cohort == "ancestral" ~ "ancestral",
                                  T ~ "descendant")) %>%
  mutate(combo = paste(substr(provenance, 1, 1), substr(cohort, 1, 1), sep = "")) %>% 
  mutate(combo = factor(combo, levels = c("ca", "sa", "cd", "sd", "hd", "kd"))) %>% 
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
  stat_summary(fun = "mean", geom = "point", 
               size = 3, aes(color = provenance, shape = cohort),
               position=position_nudge(x = 0.3)) + 
  facet_wrap(~trait, scales = "free_y", nrow = 2) + 
  scale_color_manual(values = colors) +
  ylab("trait value") +
  xlab("") + theme_bw() -> FigS7

png(here("figs_tables", "FigureS7_allTraitsRaw.png"), height = 4, width = 10,
    res = 300, units = "in")
FigS7
dev.off()

## Figure S8: stem height differs by cohort ####
mono_traits %>% 
  filter(provenance %in% c("corn", "sellman")) %>% 
  mutate(`age cohort` = case_when(`age cohort` == "ancestral" ~ "ancestral",
                                  T ~ "descendant")) %>%
  ggplot(aes(x = provenance, y = mean_tot_height)) +
  geom_boxplot(aes(color = provenance), outlier.shape = NA) +
  geom_jitter(aes(color = provenance, shape = `age cohort`), height = 0, width = 0.1,
              size = 3, alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values = colors[c(1,4)]) +
  ylab("stem height (cm)") -> height_location

png(here("figs_tables", "FigureS8_heightLoc.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
height_location
dev.off()

## Figure S9: Blue genes parameter values ####

# Parameter estimates for bMax and root:shoot from Blue Genes 2019 experiment
blue_genes <- read_rds(here("supp_data", "blue_genes_subdata.rds"))

# Fit a parabola for the aboveground biomass data
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

# Find peak biomass given a symmetric parabola
zPeak <- (zMax_for_sim + zMin_for_sim) / 2

# Calculate the predicted biomass at that elevation (bMax)
bMax <- predict(quad_mod, newdata = data.frame(elevation = zPeak))
# Convert to g / cm2
pot_area_cm2 <- pi * 5.08^2
bMax_for_sim <- bMax / pot_area_cm2

# Convert the data to g/cm2 as well
blue_genes$agb_scam_g_cm2 <- blue_genes$agb_scam / pot_area_cm2

# Figure S9a: biomass elevation parabola
blue_genes %>% 
  ggplot(aes(x = elevation*100, y = agb_scam_g_cm2)) + 
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, fullrange = T, size = 2, color = "gray47") +
  geom_hline(aes(yintercept = bMax_for_sim), linetype = "dashed", col = "purple", size = 1.5) +
  geom_vline(aes(xintercept = zMin_for_sim*100), linetype = "dotted", col = "darkgreen", size = 1.5)+
  geom_vline(aes(xintercept = zMax_for_sim*100), linetype = "dotted", col = "darkgreen", size = 1.5) +
  ylim(0, 0.22) +
  ylab(expression(paste("aboveground biomass (g/", cm^2, ")"))) +
  xlab("elevation (cm NAVD88)") +
  geom_point(aes(x = (zMin_for_sim*100 + zMax_for_sim*100)/2, y = bMax_for_sim), size = 5, color = "purple")+
  geom_point(aes(x = zMin_for_sim*100, y = 0), size = 5, color = "darkgreen") + 
  geom_point(aes(x = zMax_for_sim*100, y = 0), size = 5, color = "darkgreen") +
  theme_bw() -> fig_S9a

# Figure S9b: root-shoot
blue_genes %>%
  mutate(rs = total_bg / agb_scam) %>% 
  filter(rs < 6) %>% 
  ggplot(aes(x = rs)) +
  geom_histogram(binwidth = 0.5) +
  xlab("root-to-shoot ratio") +
  geom_vline(aes(xintercept = mean(total_bg/agb_scam)), linetype = "dashed", color = "orange", size = 1.5) +
  geom_point(aes(x = mean(total_bg/agb_scam), y = 0), size = 5, color = "orange") +
  theme_bw() + ylab("count") -> fig_S9b

png(here("figs_tables", "FigS9_BlueGenesParams.png"), height = 3.5, width = 6.5, res = 300, units = "in")
plot_grid(fig_S9a, fig_S9b, labels = "auto",
          rel_widths = c(3,2))
dev.off()
## Figure S10: trait space for MEM simulations ####

for_MEM_full <- readRDS(here("outputs/CMEM_runs", "traits_for_MEM_simulations.rds"))

# Make correlation plots between traits
# agb vs rs
cor_agb_rs <- cor(for_MEM_full$`aboveground biomass (g)`, for_MEM_full$`root:shoot ratio`)

for_MEM_full %>% 
  ggplot(aes(x = `aboveground biomass (g)`, y = `root:shoot ratio`)) +
  geom_point(alpha = 0.2, size = 2) +
  geom_smooth(method = "lm", se = F, color = "purple") +
  theme_bw(base_size = 14) +
  xlab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) -> agb_rs

# depth vs agb
cor_depth_agb <- cor(for_MEM_full$`maximum rooting depth (cm)`, for_MEM_full$`aboveground biomass (g)`)

for_MEM_full %>% 
  ggplot(aes(x = `maximum rooting depth (cm)`, y = `aboveground biomass (g)`)) +
  geom_point(alpha = 0.2, size = 2) +
  geom_smooth(method = "lm", se = F, color = "purple") +
  theme_bw(base_size = 14) +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) -> depth_agb

# rs vs depth
cor_rs_depth <- cor(for_MEM_full$`root:shoot ratio`, for_MEM_full$`maximum rooting depth (cm)`)

for_MEM_full %>% 
  ggplot(aes(x = `root:shoot ratio`, y = `maximum rooting depth (cm)`)) +
  geom_point(alpha = 0.2, size = 2) +
  geom_smooth(method = "lm", se = F, color = "purple") +
  theme_bw(base_size = 14) -> rs_depth

a10 <- ggMarginal(agb_rs, fill = "gray", color = "gray27")
b10 <- ggMarginal(depth_agb, fill = "gray", color = "gray27")
c10 <- ggMarginal(rs_depth, fill = "gray", color = "gray27")

png(here("figs_tables","FigureS10_randomdrawsMEM.png"), width = 9.22, height = 3.05,
    res = 300, units = "in")
cowplot::plot_grid(a10, b10, c10, nrow = 1)
dev.off()

