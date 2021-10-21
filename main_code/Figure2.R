# Figure 3 - Belowground Evolution

# Figure caption: 

## Preliminaries ####

# Load libraries
library(tidyverse); library(patchwork); library(here)

# Read in trait data
all_traits <- read_csv(here("data", "CompiledTraitData.csv"))
# Make a data frame of just polyculture data
poly_traits <- all_traits %>% filter(diversity == "poly")

# Read in additive expectations data from monoculture vs polyculture analysis
agb_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_agb_add.rds"))
bgb_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_bgb_add.rds"))
height_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_height_add.rds"))
width_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_width_add.rds"))
rs_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_rs_add.rds"))
density_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_density_add.rds"))
beta_additive <- readRDS(here("outputs/monoculture_polyculture", "monopoly_beta_add.rds"))

# Read in average differences by age
diffs_by_age <- readRDS(here("outputs/monoculture_polyculture", "diffs_by_age.rds"))

## Figure 2a - non-additivity across traits ####
get_avg_difference <- function(trait, additive_samples){
  observed <- pull(poly_traits[,trait])
  
  # Create matrix to hold differences between observed and predicted for each pot
  # at each iteration
  difference <- matrix(0, nrow = nrow(additive_samples$MonoPredict), ncol = ncol(additive_samples$MonoPredict))
  
  if(trait == "density"){
    for (i in 1:nrow(additive_samples$MonoPredict)){
      difference[i,] <- observed - exp(additive_samples$MonoPredict[i,])
    }
  }else{
    for (i in 1:nrow(additive_samples$MonoPredict)){
      difference[i,] <- observed - additive_samples$MonoPredict[i,]
    }
  }
  
  # Get average scaled difference across pots for each iteration
  avg_difference <- data.frame(x = rowMeans(difference, na.rm = T))/mean(observed, na.rm = T)
  return(avg_difference)
}

# Create a data frame of differences for all of the traits
diff_df <- cbind(get_avg_difference("agb", additive_samples = agb_additive),
                 get_avg_difference("bgb", additive_samples = bgb_additive),
                 get_avg_difference("rs", additive_samples = rs_additive),
                 get_avg_difference("mean_tot_height", additive_samples = height_additive),
                 get_avg_difference("mean_mid_width", additive_samples = width_additive),
                 get_avg_difference("density", additive_samples = density_additive),
                 get_avg_difference("beta", additive_samples = beta_additive))

# Add column names
colnames(diff_df) <- c("abg", "bg", "rs", "height", "width", "density", "beta")

# Create figure
tibble(diff_df) %>% 
  gather(key = trait, value = difference, abg:beta) %>% 
  mutate(trait = fct_relevel(trait, 
                             "abg", "density", "height", 
                             "width", "bg", "rs", "beta")) %>%
  mutate(color_code = case_when(trait == "rs" ~ "a",
                                trait == "beta" ~ "b",
                                T ~ "c")) %>% 
  ggplot(aes(x = trait, y = difference, color = color_code)) +
  geom_pointrange(mapping = aes(x = trait, y = difference),
                  stat = "summary",
                  fun.min = function(z) {quantile(z,0.025)},
                  fun.max = function(z) {quantile(z,0.975)},
                  fun = mean,
                  size = 0.5, shape = 1) +
  ylab("scaled difference") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray47") +
  scale_x_discrete(labels = c("ab biomass", "stem density", "stem height", "stem width",
                              "bg biomass", "root:shoot", "root parameter")) +
  annotate("text", x = 1.7, y = 0.2, label = "observed > predicted", size = 4, color = "black") +
  annotate("text", x = 1.7, y = -0.2, label = "predicted > observed", size =4, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(-0.25,0.25)) +
  theme_bw() +
  scale_color_manual(values = c("#fb9a99","#e31a1c", "gray47")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) -> Fig2a 

## Figure 2b - average non-additivity across age groups for r:s ####

diffs_by_age %>% 
  select(age, `root:shoot ratio`) %>% 
  ggplot(aes(x = age, y = `root:shoot ratio`)) +
  geom_boxplot(outlier.shape = NA, color = "#fb9a99") +
  geom_jitter(aes(shape = age), width = 0.2, alpha = 0.4, color = "#fb9a99") +
  scale_shape_manual(values = c(16,8,17)) +
  xlab("age cohort") +
  ylab("root:shoot (scaled diff)") +
  theme_bw() +
  theme(legend.position = "none") -> Fig2b

## Figure 2c - average non-additivity across age groups for root distribution ####
diffs_by_age %>% 
  select(age, `root distribution parameter`) %>% 
  ggplot(aes(x = age, y = `root distribution parameter`)) +
  geom_boxplot(outlier.shape = NA, color = "#e31a1c") +
  geom_jitter(aes(shape = age), width = 0.2, alpha = 0.4, color = "#e31a1c") +
  scale_shape_manual(values = c(16,8,17)) +
  xlab("age cohort") + 
  ylab("root parameter (scaled diff)") +
  theme_bw() +
  theme(legend.position = "none") -> Fig2c

## Join plots together ####
Fig2bc <- Fig2b + Fig2c

png(here("figs_tables", "Figure2.png"), height = 6.4, width = 6.40, units = "in", res = 300)  
Fig2a / Fig2bc + plot_layout(heights = c(3,2)) + plot_annotation(tag_levels = "a")
dev.off()  

