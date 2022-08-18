# Creates all supplemental figures for Belowground Evolution paper

## Preliminaries ####

# Load libraries
library(patchwork);library(raster); library(maps); library(cowplot);
library(ggsn); library(ggmap); library(tidyverse);
library(ggrepel); library(here); library(mvtnorm); library(GGally);
library(ggmcmc); library(ggExtra); library(rCMEM); library(gdsfmt);
library(SNPRelate)

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
blue_genes <- read_csv("supp_data/bg_data.csv")

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
  mutate(pred_year = trunc(2016 - median)) -> TableS1

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

## Figure S3: genotype PCA ####

# Read in genetic data and create usable object
#vcf_fn <- here("supp_data","94all.vcf")
#snpgdsVCF2GDS(vcf_fn, here("supp_data", "94genotypes.gds"), method = "biallelic.only")
genofile <- snpgdsOpen(here("supp_data", "94genotypes.gds"), allow.duplicate = T)
# Read in genetic data codes 
gen_data <- read_csv(here("supp_data", "genotypes_for_PCA.csv"))

# Create PCA 
pca_data <- snpgdsPCA(genofile, autosome.only = F)

# Create a tibble of the top two eigenvectors and sample ids
tibble(id = pca_data$sample.id,
       eig1 = pca_data$eigenvect[,1],
       eig2 = pca_data$eigenvect[,2],
       eig3 = pca_data$eigenvect[,3]) -> all_genotypes_pca_data

# Format genotype info data frame
names(gen_data) <- tolower(names(gen_data))
gen_data %>%
  dplyr::select(id, provenance = site, depth, cohort, expt) %>% 
  mutate(expt = ifelse(expt == 0, "no", "yes")) %>%
  mutate(expt = factor(expt, c("yes", "no"))) %>% 
  mutate(cohort = case_when(cohort == "extant" ~ "descendant (0-5cm)",
                            cohort == "descendant" ~ "descendant (0-5cm)",
                            T ~ "ancestral (10-21cm)")) %>% 
  mutate(cohort = factor(cohort, c("descendant (0-5cm)", "ancestral (10-21cm)")))-> gen_data

# Merge gen_data with SNP info
merge(gen_data, all_genotypes_pca_data) -> gen_merged

# Make a graph of the data
gen_merged %>% 
  ggplot(aes(x = eig1, y = eig2)) +
  geom_point(aes(color = expt, fill = provenance, shape = cohort),
             size = 3, stroke = 1.5, alpha = 0.7) + 
  scale_shape_manual(values = c(21, 24)) +
  theme_bw(base_size = 14) +
  xlab("axis 1 (10.9%)") +
  ylab("axis 2 (9.0%)") +
  labs(shape = "cohort",
       color = "in experiment?") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  scale_color_manual(values = c("gold", "gray47")) +
  guides(fill = guide_legend(override.aes = list(color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"))),
         color = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.position = "left",  legend.box="vertical",
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank())+
  ggtitle("All provenances")-> all_plot

# Now repeat for just the Corn genotypes alone
#vcf_corn <- here("supp_data", "Corn_only.vcf")

# Read in data and create SNP data object
#snpgdsVCF2GDS(vcf_corn, here("supp_data", "Corngenotypes.vcf"), method = "biallelic.only")
genofile_corn <- snpgdsOpen(here("supp_data", "Corngenotypes.vcf"), allow.duplicate = T)

# Create PCA
pca_data_corn <- snpgdsPCA(genofile_corn, autosome.only = F)

# Create a tibble of the top two eigenvectors and sample ids
tibble(id = pca_data_corn$sample.id,
       eig1 = pca_data_corn$eigenvect[,1],
       eig2 = pca_data_corn$eigenvect[,2]) -> corn_pca

# Merge gen_data with SNP info
merge(gen_data, corn_pca) -> gen_merged_corn

# Make a graph of the data
gen_merged_corn %>% 
  ggplot(aes(x = eig1, y = eig2)) +
  geom_point(aes(color = expt, shape = cohort),
             size = 2, stroke = 1, alpha = 0.7, fill = "#1b9e77") + 
  scale_shape_manual(values = c(21, 24)) +
  theme_bw(base_size = 14) +
  xlab("axis 1 (14.1%)") +
  ylab("axis 2 (8.4%)") +
  labs(shape = "cohort",
       color = "in experiment?") +
  scale_color_manual(values = c("gold", "gray47")) +
  guides(color = guide_legend(override.aes = list(shape = 21)))+
  theme(legend.position = "none") +
  ggtitle("Corn only")-> corn_plot

# Repeat for Sellman genotypes only
#vcf_sellman <- here("supp_data","Sellman_only.vcf")

# Read in data and create SNP data object
#snpgdsVCF2GDS(vcf_sellman, here("supp_data", "Sellmangenotypes.gds"), method = "biallelic.only")
genofile_sellman <- snpgdsOpen(here("supp_data", "Sellmangenotypes.gds"), allow.duplicate = T)

# Create PCA
pca_data_sellman <- snpgdsPCA(genofile_sellman, autosome.only = F)

# Create a tibble of the top two eigenvectors and sample ids
tibble(id = pca_data_sellman$sample.id,
       eig1 = pca_data_sellman$eigenvect[,1],
       eig2 = pca_data_sellman$eigenvect[,2]) -> sellman_pca

# Merge gen_data with SNP info
merge(gen_data, sellman_pca) -> gen_merged_sellman

# Make a graph of the data
gen_merged_sellman %>% 
  ggplot(aes(x = eig1, y = eig2)) +
  geom_point(aes(color = expt, shape = cohort),
             size = 2, stroke = 1, alpha = 0.7, fill = "#e7298a") + 
  scale_shape_manual(values = c(21, 24)) +
  theme_bw(base_size = 14) +
  xlab("axis 1 (25.7%)") +
  ylab("axis 2 (11.8%)") +
  labs(shape = "cohort",
       color = "in experiment?") +
  scale_color_manual(values = c("gold", "gray47")) +
  guides(color = guide_legend(override.aes = list(shape = 21)))+
  theme(legend.position = "none") +
  ggtitle("Sellman only")-> sellman_plot

# Bring plots together
design <- c(
  patchwork::area(1, 1, 2, 2),
  patchwork::area(1, 3, 1, 3),
  patchwork::area(2, 3, 2, 3)
)

png("figs_tables/FigS3_SNP_PCA.png", height = 7, width = 11, res = 300, units = "in")
all_plot + corn_plot + sellman_plot +
  plot_layout(guides = "collect", design = design) +
  plot_annotation(tag_levels = "a")
dev.off()

## Figure S4: random intercepts of monoculture models ####
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

png(here("figs_tables", "FigS4_genotype.png"), height = 8.6, width = 10, res = 300, units = "in")
left_panel + right_panel
dev.off() 

## Figure S5: all traits by cohort with raw data ####
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
  xlab("") + theme_bw() -> FigS5

png(here("figs_tables", "FigS5_allTraitsRaw.png"), height = 4, width = 10,
    res = 300, units = "in")
FigS5
dev.off()

## Figure S6: root-to-shoot ratio differs by provenance and cohort ####
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

png(here("figs_tables", "FigS6_rsLocAge.png"), height = 2.7, width = 5.4,
    res = 300, units = "in")
rs_location + rs_age + plot_layout(guides = "collect")
dev.off()

## Figure S7: stem width differs by cohort ####
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

png(here("figs_tables", "FigS7_widthAge.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
width_age
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

png(here("figs_tables", "FigS8_heightLoc.png"), height = 2.7, width = 3.5,
    res = 300, units = "in")
height_location
dev.off()

## Figure S9: monopoly diffs by age ####
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
        legend.justification = c(1.6, 0)) -> fig_S9

png(here("figs_tables","FigS9_monopoly_byAge.png"), height = 4, width = 10, units = "in", res = 300)
fig_S9
dev.off()

# Also test for differences across age cohorts
anova(lm(`aboveground biomass (g)` ~ age, data = diffs_by_age)) # ns
anova(lm(`stem density` ~ age, data = diffs_by_age)) # ns
anova(lm(`mean stem height (cm)` ~ age, data = diffs_by_age)) # ns
anova(lm(`mean stem width (mm)` ~ age, data = diffs_by_age)) # ns
anova(lm(`belowground biomass (g)` ~ age, data = diffs_by_age)) # ns
anova(lm(`root:shoot ratio` ~ age, data = diffs_by_age)) # ns
anova(lm(`root distribution parameter` ~ age, data = diffs_by_age)) # .

## Figure S10: Expt 2 parameter values ####

# Subset blue genes data for root-to-shoot ratio and aboveground biomass
# parameter estimation -- do not include genotypes from Blackwater or
# competition pots that included Spartina patens
blue_genes %>% 
  filter(provenance != "blackwater" & comp == 0 & co2 == "ambient") -> blue_genes_sub

# Convert to g / cm2
pot_area_cm2 <- pi * 5.08^2
blue_genes_sub$agb_scam_g_cm2 <- blue_genes_sub$agb_scam / pot_area_cm2

# Read in parameter data for simulations
bg_params <- read_rds(here("outputs", "bg_params_plotting.rds"))

# Figure S10a: biomass elevation parabola
blue_genes_sub %>% 
  ggplot(aes(x = elevation*100, y = agb_scam_g_cm2)) + 
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, fullrange = T, size = 2, color = "gray47") +
  geom_hline(aes(yintercept = bg_params$bMax), linetype = "dashed", col = "purple", size = 1.5) +
  geom_vline(aes(xintercept = bg_params$zMin*100), linetype = "dotted", col = "darkgreen", size = 1.5)+
  geom_vline(aes(xintercept = bg_params$zMax*100), linetype = "dotted", col = "darkgreen", size = 1.5) +
  ylim(0, 0.18) +
  ylab(expression(paste("aboveground biomass (g/", cm^2, ")"))) +
  xlab("elevation (cm NAVD88)") +
  geom_point(aes(x = (bg_params$zMin*100 + bg_params$zMax*100)/2, y = bg_params$bMax), size = 5, color = "purple")+
  geom_point(aes(x = bg_params$zMin*100, y = 0), size = 5, color = "darkgreen") + 
  geom_point(aes(x = bg_params$zMax*100, y = 0), size = 5, color = "darkgreen") +
  theme_bw(base_size = 14) -> fig_S10a 

# Figure S10b: root-shoot

blue_genes_sub %>%
  mutate(rs = total_bg / agb_scam) %>% 
  filter(rs < 5) %>% 
  ggplot(aes(x = rs)) +
  geom_histogram(binwidth = 0.5) +
  xlab("root-to-shoot ratio") +
  geom_vline(aes(xintercept = mean(total_bg/agb_scam)), linetype = "dashed", color = "orange", size = 1.5) +
  geom_point(aes(x = mean(total_bg/agb_scam), y = 0), size = 5, color = "orange") +
  theme_bw(base_size = 14) + ylab("count") ->fig_S10b

png(here("figs_tables", "FigS10_BlueGenesParams.png"), height = 3.5, width = 6.5, res = 300, units = "in")
plot_grid(fig_S10a, fig_S10b, labels = "auto",
          rel_widths = c(3,2))
dev.off()
## Figure S11: trait space for CMEM simulations ####

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

a11 <- ggMarginal(agb_rs, fill = "gray", color = "gray27")
b11 <- ggMarginal(depth_agb, fill = "gray", color = "gray27")
c11 <- ggMarginal(rs_depth, fill = "gray", color = "gray27")

png(here("figs_tables","FigS11_randomdrawsMEM.png"), width = 9.22, height = 3.05,
    res = 300, units = "in")
cowplot::plot_grid(a11, b11, c11, nrow = 1)
dev.off()


## Figure S12: Expt 2 analysis: genotype-level parabolas ####
# Filter data to Corn Island only for genotypes that have at least 20 reps
blue_genes %>% 
  filter(provenance == "corn") %>% 
  group_by(genotype) %>%
  filter(n() >= 20)-> quad_test_data

# Create a vector of genotype names for loop procedure later
sort(unique(quad_test_data$genotype)) -> quad_genotypes

# Create storage to hold coefficients from quadratic regressions of each
# genotype
param_storage <- matrix(NA, nrow = length(quad_genotypes), ncol = 3)

# Fit quadratic regression to agb data for each genotype and store coefficients
for (i in 1:length(quad_genotypes)){
  temp_data <- quad_test_data %>% filter(genotype == quad_genotypes[i])
  temp_mod <- lm(agb_scam ~ elevation + I(elevation^2), data = temp_data)
  params <- as.numeric(coef(temp_mod))
  param_storage[i,] <- c(params[1], params[2], params[3])
}

# Create function to convert coefficients to x-intercepts
quadraticRoots <- function(a, b, c) {
  discriminant <- (b^2) - (4*a*c)
  x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
  x_int_neg <- (-b - sqrt(discriminant)) / (2*a)
  xints <- c(x_int_plus, x_int_neg)
  return(xints)
}

# Calculate x-intercepts for all genotypes
matrix(quadraticRoots(param_storage[,3], param_storage[,2], param_storage[,1]),
       nrow = length(quad_genotypes), ncol = 2)-> x_ints

# Find elevation where peak biomass is given a symmetric parabola
zPeak <- (x_ints[,1] + x_ints[,2]) / 2

# Calculate the predicted biomass at that elevation for each genotype (bMax)
bMax <- NULL
for (i in 1:length(quad_genotypes)){
  temp_data <- quad_test_data %>% filter(genotype == quad_genotypes[i])
  temp_mod <- lm(agb_scam ~ elevation + I(elevation^2), data = temp_data)
  bMax[i] <- predict(temp_mod, newdata = data.frame(elevation = zPeak[i]))
}

# Collect number of reps per genotype
quad_test_data %>% 
  group_by(genotype) %>% 
  summarize(count = n()) %>% 
  pull(count) -> reps

# Put all parameter info together in a table
tibble(lower_x = x_ints[,1],
       upper_x = x_ints[,2],
       b_max = bMax,
       n = reps,
       genotype = quad_genotypes,
       range = x_ints[,2] - x_ints[,1]) %>% 
  # Create age cohort column
  mutate(age = ifelse(substr(genotype, 2, 2) == "a",
                      "ancestral", "modern"))-> quad_params_by_genotype

# Get summary stats by age cohort for each parameter  
quad_params_by_genotype %>% 
  group_by(age) %>% 
  summarize(range = mean(range),
            bmax = mean(b_max), 
            lower = mean(lower_x),
            upper = mean(upper_x),
            n = n()) %>% 
  # Convert biomass to g/cm2 and elevations to cm
  mutate(bmax_gcm2 = bmax / (pi*5.08^2),
         lower_cm = lower*100,
         upper_cm = upper*100,
         range_cm = range*100,
         peak_cm = lower_cm + (upper - lower)*100/2) %>% 
  dplyr::select(age, n, bmax_gcm2, lower_cm, upper_cm, range_cm, peak_cm) -> summary_params

# Create a plot of quadratic curves for each genotype (color by age cohort)
quad_test_data %>%
  mutate(age = ifelse(age == "ancestral", "ancestral", "descendant")) %>% 
  mutate(agb = agb_scam / (pi*5.08^2)) %>% 
  ggplot(aes(x = elevation*100, y = agb, group = genotype, color = age)) +
  geom_point(alpha = 0.5) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x + I(x^2),
            alpha = 0.3, se = F, fullrange = T, size = 1.1) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  xlim(-10, 70) +
  ylim(0, 0.16) +
  theme(legend.position=c(0.8,0.9),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(face = 2)) -> genotype_plot

png(here("figs_tables", "FigS12_Expt2Parabolas.png"), width = 4.7, height = 4.2,
    res = 300, units = "in")
genotype_plot
dev.off()

## Figure S13: Expt 2 analysis: parameter scenarios ####
# Get tidal data
tides <- read_csv("~/Git/belowground_evolution/supp_data/tides_2018.csv")
mean(tides$`MSL (m)`) -> msl
mean(tides$`MHW (m)`) -> mhw

# Set up storage
run_store <- matrix(NA, nrow = 6, ncol = 100)
carbon_store <- matrix(NA, nrow = 6, ncol = 100)

# Set up parameter scenarios
# Scenario 1 - ancestral bMax, average values for z_min and z_max
tibble(scenario = c(1,2,3,4,5,6),
       bmax = c(summary_params$bmax_gcm2[1], summary_params$bmax_gcm2[2], mean(summary_params$bmax_gcm2),
                mean(summary_params$bmax_gcm2),summary_params$bmax_gcm2[1], summary_params$bmax_gcm2[2]),
       zmin = c(mean(summary_params$lower_cm), mean(summary_params$lower_cm), summary_params$lower_cm[1],
                summary_params$lower_cm[2], summary_params$lower_cm[1], summary_params$lower_cm[2]),
       zmax = c(mean(summary_params$upper_cm), mean(summary_params$upper_cm), summary_params$upper_cm[1],
                summary_params$upper_cm[2], summary_params$upper_cm[1], summary_params$upper_cm[2])) -> model_scenarios

for (i in 1:6){
  mem_out <- rCMEM::runCohortMem(startYear=2020, relSeaLevelRiseInit=0.34, relSeaLevelRiseTotal=34,
                                 initElv=22.6, meanSeaLevel=msl*100,
                                 meanHighWaterDatum=mhw*100, suspendedSediment=3e-05,
                                 lunarNodalAmp=0, bMax = model_scenarios$bmax[i], 
                                 zVegMin = model_scenarios$zmin[i], zVegMax=model_scenarios$zmax[i], zVegPeak=NA,
                                 plantElevationType="orthometric", rootToShoot = 1.8,
                                 rootTurnover=0.55, rootDepthMax=30, omDecayRate=0.8,
                                 recalcitrantFrac=0.2, captureRate = 2.8)
  run_store[i,] <- mem_out$annualTimeSteps$surfaceElevation
  
  mem_out$cohorts %>% 
    mutate(loi = (fast_OM + slow_OM + root_mass) / (fast_OM + slow_OM + root_mass + mineral),
           perc_C = 0.4*loi + 0.0025*loi^2,
           layer_C = (fast_OM + slow_OM + root_mass)*perc_C) %>% 
    group_by(year) %>% 
    summarize(total_C = sum(layer_C)) %>% pull(total_C) -> carbon_store[i,]
  print(i)
}

# Calculate average accretion rates
init_elev <- 22.6
avg_accretion_rate <- (run_store[,80] - init_elev) / 80

# Calculate average carbon accumulation rates
avg_C_accum_rate <- (carbon_store[,80] - carbon_store[,1]) / 80

# Create a data frame to hold all of that information
tibble(scenario = c("ancestral (bMax)", "descendant (bMax)",
                    "ancestral (range)", "descendant (range)",
                    "ancestral (bMax + range)", "descendant (bMax + range)"),
       `vert. accretion rate` = avg_accretion_rate * 10,
       # unit conversion for carbon accumulation 
       `C accum. rate` = avg_C_accum_rate * 1e-6 / 1e-8) 

# Create second summary plot of mean predictions and parameter values for each
# age cohort using the 3 scenarios outlined about

# Convert Z_min, Z_max, and b_max back to a,b,c for ancestral and modern
quadraticABC <- function(zVegMin, bMax, zVegMax){
  zVegPeak <- (zVegMax - zVegMin)/2
  a <- -((-zVegMin * bMax - zVegMax * bMax) / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  b <- -(bMax / ((zVegMin - zVegPeak) * (-zVegMax + zVegPeak)))
  c <- (zVegMin * zVegMax * bMax) / ((zVegMin - zVegPeak) * (zVegMax - zVegPeak))
  return(c(a,b,c))
}

# Scenario 1 - vary bmax but don't vary intercepts
params_out1 <- quadraticABC(rep(mean(summary_params$lower_cm),2), summary_params$bmax_gcm2, rep(mean(summary_params$upper_cm),2))
x_pred <- seq(-10, 70, 0.1)
y_pred_anc1 <- x_pred * params_out1[1] + x_pred^2 * params_out1[3] + params_out1[5]
y_pred_mod1 <- x_pred * params_out1[2] + x_pred^2 * params_out1[4] + params_out1[6] 

tibble(elevation = rep(x_pred,2),
       agb = c(y_pred_anc1, y_pred_mod1),
       age = rep(c("ancestral", "descendant"), each = length(x_pred))) %>% 
  ggplot(aes(x = elevation, y = agb, color = age)) +
  geom_line(size = 2) +
  ylim(0, 0.08) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  # Add Z_min for each age cohort
  geom_point(aes(x = mean(summary_params$lower_cm), y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = mean(summary_params$lower_cm), y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = mean(summary_params$upper_cm), y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = mean(summary_params$upper_cm), y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = mean(summary_params$peak_cm), y = summary_params$bmax_gcm2[1]), size = 5, color = "#D55E00") +
  geom_point(aes(x = mean(summary_params$peak_cm), y = summary_params$bmax_gcm2[2]), size = 5, color = "#0072B2") +
  theme_bw(base_size = 16) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.position=c(0.8,0.9),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.text = element_text(face = 2))+
  ggtitle("Scenario 1: vary bMax only") -> cohort_plot1

# Scenario 2 - vary intercepts but don't vary bmax
params_out2 <- quadraticABC(summary_params$lower_cm, rep(mean(summary_params$bmax_gcm2),2), summary_params$upper_cm)
y_pred_anc2 <- x_pred * params_out2[1] + x_pred^2 * params_out2[3] + params_out2[5]
y_pred_mod2 <- x_pred * params_out2[2] + x_pred^2 * params_out2[4] + params_out2[6] 

tibble(elevation = rep(x_pred,2),
       agb = c(y_pred_anc2, y_pred_mod2),
       age = rep(c("ancestral", "modern"), each = length(x_pred))) %>% 
  ggplot(aes(x = elevation, y = agb, color = age)) +
  geom_line(size = 2) +
  ylim(0, 0.08) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  # Add Z_min for each age cohort
  geom_point(aes(x = summary_params$lower_cm[1], y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$lower_cm[2], y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = summary_params$upper_cm[1], y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$upper_cm[2], y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = summary_params$peak_cm[1], y = mean(summary_params$bmax_gcm2)), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$peak_cm[2], y = mean(summary_params$bmax_gcm2)), size = 5, color = "#0072B2") +
  theme_bw(base_size = 16) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme(legend.position = "none", plot.title = element_text(size = 18, face = "bold")) +
  ggtitle("Scenario 2: vary zMin and zMax only") -> cohort_plot2

# Scenario 3 - vary intercepts AND bmax
params_out3 <- quadraticABC(summary_params$lower_cm, summary_params$bmax_gcm2, summary_params$upper_cm)
y_pred_anc3 <- x_pred * params_out3[1] + x_pred^2 * params_out3[3] + params_out3[5]
y_pred_mod3 <- x_pred * params_out3[2] + x_pred^2 * params_out3[4] + params_out3[6] 

tibble(elevation = rep(x_pred,2),
       agb = c(y_pred_anc3, y_pred_mod3),
       age = rep(c("ancestral", "modern"), each = length(x_pred))) %>% 
  ggplot(aes(x = elevation, y = agb, color = age)) +
  geom_line(size = 2) +
  ylim(0, 0.08) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  # Add Z_min for each age cohort
  geom_point(aes(x = summary_params$lower_cm[1], y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$lower_cm[2], y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = summary_params$upper_cm[1], y = 0), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$upper_cm[2], y = 0), size = 5, color = "#0072B2") +
  geom_point(aes(x = summary_params$peak_cm[1], y = summary_params$bmax_gcm2[1]), size = 5, color = "#D55E00") +
  geom_point(aes(x = summary_params$peak_cm[2], y = summary_params$bmax_gcm2[2]), size = 5, color = "#0072B2") +
  theme_bw(base_size = 16) +
  xlab("elevation (cm NAVD88)") +
  ylab(expression(paste("aboveground biomass (g ", cm^-2, ")"))) +
  theme(legend.position = "none", plot.title = element_text(size = 18, face = "bold")) +
  ggtitle("Scenario 3: vary bMax, zMin, and zMax") -> cohort_plot3

png(here("figs_tables", "FigS13_Expt2Scenarios.png"), height = 12, width = 6, res = 300, units = "in")
cohort_plot1 / cohort_plot2 / cohort_plot3 
dev.off()
