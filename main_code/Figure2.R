# Belowground evolution Figure 2

## Preliminaries ####

# Load libraries
library(tidyverse); library(cowplot); library(here); library(ggmcmc); library(lme4)

# Set colors for plotting
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Read in Bayesian regression model output from corn-sellman models
cs_agb_out <- read_rds(here("outputs/corn_sellman_models", "agb_csmodel.rds"))
cs_bgb_out <- read_rds(here("outputs/corn_sellman_models", "bgb_csmodel.rds"))
cs_width_out <- read_rds(here("outputs/corn_sellman_models", "width_csmodel.rds"))
cs_height_out <- read_rds(here("outputs/corn_sellman_models", "height_csmodel.rds"))
cs_density_out <- read_rds(here("outputs/corn_sellman_models", "density_csmodel.rds"))
cs_rs_out <- read_rds(here("outputs/corn_sellman_models", "rs_csmodel.rds"))
cs_beta_out <- read_rds(here("outputs/corn_sellman_models", "beta_csmodel.rds"))

# Read in raw trait data
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))

# Filter to get just corn and sellman
all_traits %>% 
  filter(location %in% c("corn", "sellman")) -> cs_traits

# Create a model template
model_template_cs <- lmer(agb ~ ln_depth + ic_weight + frame +
                            location + age + (1|genotype), data = cs_traits)


## Figure 2a - proportional variance explained by terms in linear model ####

# Create a function to calculate the proportion of variance that comes from each
# of the different components (experimental conditions, provenance, age,
# genotype, residual)

calculate_prop_var <- function(coda_object, trait, model_template){
  
  # First calculate variation due to genotype and residual (Nakagawa et al. for
  # Poisson GLMM)
  if(trait == "density"){
    ggs(coda_object) %>% 
      spread(key = Parameter, value = value) %>%
      mutate(sigma2.res = log(1 + 1/mean(cs_traits$density))) %>% 
      mutate(sigma2.int = sigma.int^2) %>% 
      dplyr::select(sigma2.res, sigma2.int) -> sigma_out
  }else{
    ggs(coda_object) %>% 
      filter(Parameter %in% c("sigma.int", "sigma.res")) %>% 
      spread(key = Parameter, value = value) %>% 
      mutate(sigma2.int = sigma.int^2,
             sigma2.res = sigma.res^2) %>% 
      dplyr::select(sigma2.int, sigma2.res) -> sigma_out
  }
  
  # Pull regression coefficient for provenance
  ggs(coda_object) %>% 
    filter(Parameter == c("beta[6]")) %>% 
    pull(value) -> beta_6
  
  # Pull regression coefficient for age cohort
    ggs(coda_object) %>% 
      filter(Parameter == c("beta[7]")) %>% 
      pull(value) -> beta_7
  
  # Pull regression coefficients for initial weight, peat depth, and frame
  ggs(coda_object) %>% 
    spread(key = Parameter, value = value) %>% 
    dplyr::select(`beta[1]`, `beta[2]`,`beta[3]`,`beta[4]`,`beta[5]`) -> beta_15
  
  beta_15 <- as.matrix(as.data.frame(beta_15))
  
  # Create vectors of binary indicators for provenance and age
  loc_ind <- rep(0,42)
  loc_ind[which(cs_traits$location == "corn")] <- 1
  
  age_ind <- rep(0,42)
  age_ind[which(cs_traits$age == "ancestral")] <- 1
  
  # Predict variance due to location and age (note this does not account for
  # covariances between location and age and initial conditions so it is an
  # approximation)
  pred6 <- NULL
  pred7 <- NULL
  for (i in 1:length(beta_6)){
    pred6[i] <- var(beta_6[i] * loc_ind)
    pred7[i] <- var(beta_7[i] * age_ind)
    }
  
  # Predict variance due to experimental factors (initial conditions)
  pred_15 <- matrix(NA, nrow = nrow(beta_15), ncol = 42)
  X <- as.matrix(model.matrix(model_template))[1:42, 2:6]
  for (i in 1:nrow(beta_15)){
    pred_15[i,] <- var(X %*% beta_15[i,])
  }
  
  # Calculate proportion of variance due to provenance, location, genotype, residual, and initial conditions 
  fixed6 <- mean(pred6 / (sigma_out$sigma2.int + pred6 + pred7 + sigma_out$sigma2.res + pred_15))
  fixed7 <- mean(pred7 / (sigma_out$sigma2.int + pred6 + pred7 + sigma_out$sigma2.res + pred_15))
  genotype <- mean(sigma_out$sigma2.int / (sigma_out$sigma2.int + pred6 + pred7 + sigma_out$sigma2.res + pred_15)) 
  resid <- mean(sigma_out$sigma2.res / (sigma_out$sigma2.int + pred6 + pred7 + sigma_out$sigma2.res + pred_15)) 
  other <- mean(pred_15 / (sigma_out$sigma2.int + pred6 + pred7 + sigma_out$sigma2.res + pred_15)) 
  
  # Return as a list
  return(list(fixed6 = fixed6, fixed7 = fixed7, genotype = genotype, resid = resid, other = other))

}

# Set frame to be a factor
cs_traits$frame <- as.factor(cs_traits$frame)

# Create model template
model_template_cs <- lmer(agb ~ ln_depth + ic_weight + frame +
                            location + age + (1|genotype), data = cs_traits)

prop_biomass <- as.numeric(calculate_prop_var(cs_agb_out, "agb", model_template_cs))
prop_width <- as.numeric(calculate_prop_var(cs_width_out, "width", model_template_cs))
prop_bgb <- as.numeric(calculate_prop_var(cs_bgb_out, "bgb", model_template_cs))
prop_rs <- as.numeric(calculate_prop_var(cs_rs_out, "rs", model_template_cs))
prop_density <- as.numeric(calculate_prop_var(cs_density_out, "density", model_template_cs))
prop_height<- as.numeric(calculate_prop_var(cs_height_out, "height", model_template_cs))

# Calculate beta separately because no random effect of genotype

# Pull out variance due to genotype
ggs(cs_beta_out) %>% 
  filter(Parameter == "sigma.res") %>% 
  spread(key = Parameter, value = value) %>% 
  mutate(sigma2.res = sigma.res^2) %>% 
  dplyr::select(sigma2.res) -> sigma_beta

# Pull out regression coefficient for location
ggs(cs_beta_out) %>% 
  filter(Parameter == "beta[6]") %>% 
  pull(value)-> beta6_beta

# Pull out regression coefficient for provenance
ggs(cs_beta_out) %>% 
  filter(Parameter == "beta[7]") %>% 
  pull(value)-> beta7_beta

# Pull out regression coefficients due to initial conditions
ggs(cs_beta_out) %>% 
  spread(key = Parameter, value = value) %>% 
  dplyr::select(`beta[1]`, `beta[2]`,`beta[3]`,`beta[4]`,`beta[5]`) -> beta_15_beta

beta_15_beta <- as.matrix(as.data.frame(beta_15_beta))

# Calculate variance due to initial conditions
pred_15_beta <- matrix(NA, nrow = length(beta6_beta), ncol = 42)
X <- as.matrix(model.matrix(model_template_cs))[1:42, 2:6]
for (i in 1:length(beta6_beta)){
  pred_15_beta[i,] <- var(X %*% beta_15_beta[i,])
}

# Calculate variance due to provenance and age cohort
pred_beta6 <- NULL
pred_beta7 <- NULL

loc_ind <- as.numeric(as.factor(cs_traits$location)) - 1
age_ind <- as.numeric(as.factor(cs_traits$age)) - 1

for (i in 1:length(beta6_beta)){
  pred_beta6[i] <- var(beta6_beta[i] * loc_ind)
  pred_beta7[i] <- var(beta7_beta[i] * age_ind)
}

# Calculate proportion of variance due to provenance, location, residual, and initial conditions 
fixed6 <- mean(pred_beta6 / (sigma_beta$sigma2.res + pred_beta6 + pred_beta7 + pred_15_beta))
fixed7 <- mean(pred_beta7 / (sigma_beta$sigma2.res + pred_beta6 + pred_beta7 + pred_15_beta))
residual <- mean(sigma_beta$sigma2.res / (sigma_beta$sigma2.res + pred_beta6 + pred_beta7 + pred_15_beta))
other <- mean(pred_15_beta / (sigma_beta$sigma2.res + pred_beta6 + pred_beta7 + pred_15_beta))
prop_beta <- c(fixed6, fixed7, 0, residual, other)

# Now put together and graph
tibble(trait = c("abg", "density", "height", "width", "bgb", "rs", "beta"),
                     cohort = c(prop_biomass[2], prop_density[2], prop_height[2], prop_width[2],prop_bgb[2], prop_rs[2], prop_beta[2]),
                     provenance = c(prop_biomass[1], prop_density[1], prop_height[1], prop_width[1],prop_bgb[1], prop_rs[1], prop_beta[1]),
                     genotype = c(prop_biomass[3], prop_density[3], prop_height[3], prop_width[3],prop_bgb[3], prop_rs[3], prop_beta[3]),
                     residual = c(prop_biomass[4], prop_density[4], prop_height[4], prop_width[4],prop_bgb[4], prop_rs[4], prop_beta[4]),
                     other = c(prop_biomass[5], prop_density[5], prop_height[5], prop_width[5],prop_bgb[5], prop_rs[5], prop_beta[5])) %>% 
  gather(key = component, value = value, other:cohort) %>% 
  mutate(trait = fct_relevel(trait, 
                             "beta", "rs", "bgb", 
                             "width", "height", "density", "abg")) %>% 
  mutate(component = fct_relevel(component, "residual", "other", "cohort", "provenance", "genotype")) %>% 
  ggplot(aes(x = trait, y = value, fill = component)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("gray67", "gray87","#08519c","#3182bd","#9ecae1")) +
  coord_flip() +
  scale_x_discrete(labels = c("root parameter","root:shoot", "bg biomass", "stem width",
                              "stem height", "stem density", "ag biomass")) +
  ylab("proportion of variance") + theme_classic()+ theme(legend.position = "top", legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) -> stack_plot


## Figure 2b - provenance and age cohort mediate root distribution ####

# Predict cumulative root distribution for all corn-sellman monoculture pots
depth_roots <- seq(0, 50, length.out = 1000)
pred_store <- matrix(NA, nrow = length(depth_roots), ncol = nrow(cs_traits))
for(j in 1:nrow(cs_traits)){
  for (i in 1:length(depth_roots)){
    pred_store[i,j] <- 1 - cs_traits$beta[j] ^ (depth_roots[i])
  }
}

# Also get averages for each provenance by age cohort combination
cs_traits %>% 
  group_by(location, age) %>% 
  summarize(mean_beta = mean(beta)) -> beta_averages

# Calculate average rooting depth at 95% cumulative probability for each
# provenance by age cohort
depths_averages <- rep(NA, 4)

depths_averages[1]<- depth_roots[which(abs(1 - beta_averages$mean_beta[1]^depth_roots-0.95)==min(abs(1 - beta_averages$mean_beta[1]^depth_roots-0.95)))]
depths_averages[2]<- depth_roots[which(abs(1 - beta_averages$mean_beta[2]^depth_roots-0.95)==min(abs(1 - beta_averages$mean_beta[2]^depth_roots-0.95)))]
depths_averages[3]<- depth_roots[which(abs(1 - beta_averages$mean_beta[3]^depth_roots-0.95)==min(abs(1 - beta_averages$mean_beta[3]^depth_roots-0.95)))]
depths_averages[4]<- depth_roots[which(abs(1 - beta_averages$mean_beta[4]^depth_roots-0.95)==min(abs(1 - beta_averages$mean_beta[4]^depth_roots-0.95)))]

# Convert predicted probabilities to a data frame and format
as_tibble(pred_store) -> pred_store
names(pred_store) <- as.character(cs_traits$pot)
pred_store %>% 
  mutate(depth = depth_roots) %>% 
  gather(key = pot, value = cumsum, `49`:`96`) %>% 
  mutate(pot = as.numeric(pot)) -> pred_store

# Merge with trait data
merge(pred_store, cs_traits, by = "pot") -> plot_data 

# Create plot with lines for each pot colored by provenance and with linetypes
# for age cohort
plot_data %>% 
  ggplot(aes(x = depth.x, y = cumsum, group = pot, color = location, linetype = age)) + 
  geom_line(alpha = 0.2, size = 0.8) +
  geom_hline(aes(yintercept = 0.95)) +
  scale_x_reverse() +
  scale_y_reverse()+
  coord_flip() +
  ylab("cumulative probability") +
  xlab("depth below soil surface (cm)") +
  scale_color_manual(values = colors[c(1,4)]) +
  geom_point(aes(y = 0.95, x = depths_averages[1]), color = colors[1], size = 4, shape = 16) +
  geom_point(aes(y = 0.95, x = depths_averages[2]), color = colors[1], size = 4, shape = 17) +
  geom_point(aes(y = 0.95, x = depths_averages[3]), color = colors[4], size = 4, shape = 16)+
  geom_point(aes(y = 0.95, x = depths_averages[4]), color = colors[4], size = 4, shape = 17) +
  theme_bw() +
  theme(legend.position = "none") -> b_main

# Create inset graph to show differences in parameter values
cs_traits %>%
  mutate(age = case_when(age == "modern" ~ "descendant",
                         T ~ age)) %>%
  mutate(age = factor(age, levels = c("ancestral", "mix", "descendant"))) %>% 
  mutate(location = case_when(location == "corn" ~ "Corn Island",
                              T ~ "Sellman Creek")) %>% 
  ggplot(aes(x = age, y = beta)) +
  geom_boxplot(aes(color = location), outlier.shape = NA) +
  geom_jitter(aes(shape = age, color = location), height = 0, width = 0.1, alpha = 0.4, size = 2) +
  facet_wrap(~location) +
  scale_color_manual(values = colors[c(1,4)]) + theme_bw() +
  theme(legend.position = "none",
        plot.background = element_rect(colour = "gray47", fill=NA, size=1)) +
  ylab(expression(paste("root distribution parameter (", beta, ")"))) +
  xlab("age cohort") -> b_inset

ggdraw() +
  draw_plot(b_main) +
  draw_plot(b_inset, x = 0.25, y = 0.15, width = 0.7, height = 0.4) -> fig2_beta_updated


## Bring plots together ####
png(here("figs_tables", "Figure2.png"), height = 5.5, width = 11, res = 300, units = "in")
plot_grid(stack_plot, fig2_beta_updated, labels = "auto", rel_widths = c(2,2))
dev.off()

## For in-text ####
# Calculate the proportions that are due to genotype, ecotype, and age combined for each trait
list(prop_biomass,
     prop_density,
     prop_height,
     prop_width,
     prop_bgb,
     prop_rs,
     prop_beta) -> all_props

sapply(all_props, function(x) sum(x[1:3])) -> all_genotype_prop
range(all_genotype_prop)
# Not accounting for within-species variation results in missing out on
# explaining 16.5-52.6% of observed variation in traits