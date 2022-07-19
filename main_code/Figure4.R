# Belowground evolution Figure 4

## Preliminaries ####

# Load libraries
library(tidyverse); library(cowplot); library(here)

# Read in outputs from CMEM runs where agb, r:s, and rooting depth were
# manipulated
cmem_full <- read_rds(here("outputs/CMEM_runs", "CMEM_predictions_full.rds"))
# Read in outputs from CMEM runs where only agb was manipulated
cmem_agb <- read_rds(here("outputs/CMEM_runs", "CMEM_predictions_agb_only.rds"))
# Read in average accretion and C accumulation rates from full model
cmem_rates <- read_rds(here("outputs/CMEM_runs", "CMEM_rates_full.rds"))
# Read in average accretion and C accumulation rates from full model by cohort
cmem_rates_cohort <- read_rds(here("outputs/CMEM_runs", "CMEM_rates_full_cohort.rds"))

# Set colors for plotting
colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
## Figure 4a - surface elevation predictions for AG + BG simulations ####

# Collect elevation values for each age-location cohort at the last year of the
# simulation
cmem_full %>% 
  filter(iteration %in% c("corn-ancestral", "sellman-ancestral",
                          "corn-modern", "sellman-modern") & year == "2100") %>% 
  pull(value) -> end_points

cmem_full %>% 
  mutate(year = as.numeric(year)) %>%
  mutate(color_code = case_when(iteration == "corn-ancestral" ~ 1,
                                iteration == "sellman-ancestral" ~ 2,
                                iteration == "corn-modern" ~ 3,
                                iteration == "sellman-modern" ~ 4,
                                T ~ 0),
         size_code = case_when(iteration %in% c("corn-ancestral", "sellman-ancestral",
                                                "corn-modern", "sellman-modern") ~ "big",
                               T ~ "small")) %>%
  ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) +
  geom_line(aes(alpha = size_code)) +
  ylab("marsh elevation (cm NAVD88)") +
  scale_color_manual(values = c("gray11", colors[1], colors[4], colors[1], colors[4])) +
  scale_size_manual(values = c(1.5, 0.8)) +
  scale_alpha_manual(values = c(0.8, 0.1)) +
  geom_point(aes(x = 2100, y = end_points[1]), color = colors[1], size = 3) +
  geom_point(aes(x = 2100, y = end_points[2]), color = colors[4], size = 3) +
  geom_point(aes(x = 2100, y = end_points[3]), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = 2100, y = end_points[4]), color = colors[4], size = 3, shape = 17) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(22.5,35) -> Fig4_panelA

## Figure 4b - average accretion rate histogram ####
tibble(acc_rate = cmem_rates$avg_acc*10) %>%
  ggplot(aes(x = acc_rate)) +
  geom_histogram(bins=10, color = "gray27", fill = "white") +
  xlab(expression(paste("vertical accretion rate (mm ",yr^-1,")"))) +
  geom_point(aes(x = cmem_rates_cohort$acc_v[1], y = 0), color = colors[1], size = 3) +
  geom_point(aes(x = cmem_rates_cohort$acc_v[2], y = 0), color = colors[4], size = 3) +
  geom_point(aes(x = cmem_rates_cohort$acc_v[3], y = 0), color = colors[1], size = 3, shape = 17)+
  geom_point(aes(x = cmem_rates_cohort$acc_v[4], y = 0), color = colors[4], size = 3, shape = 17)+
  geom_segment(aes(x = cmem_rates_cohort$acc_v[1], y = 0, xend = cmem_rates_cohort$acc_v[1], yend = Inf),
               color = colors[1], size = 1.5, linetype = "dotted") +
  geom_segment(aes(x = cmem_rates_cohort$acc_v[2], y = 0, xend = cmem_rates_cohort$acc_v[2], yend = Inf),
               color = colors[4], size = 1.5, linetype = "dotted") + 
  geom_segment(aes(x = cmem_rates_cohort$acc_v[3], y = 0, xend = cmem_rates_cohort$acc_v[3], yend = Inf),
               color = colors[1], size = 1.5, linetype = "dotted") +
  geom_segment(aes(x = cmem_rates_cohort$acc_v[4], y = 0, xend = cmem_rates_cohort$acc_v[4], yend = Inf),
               color = colors[4], size = 1.5, linetype = "dotted") +
  theme_bw() +
  ylab("count") -> Fig4_panelB

## Figure 4c - average carbon accumulation rate histogram ####

tibble(acc_rate = cmem_rates$avg_C * 1e-6 / 1e-8) %>%
  ggplot(aes(x = acc_rate)) +
  geom_histogram(bins=10, color = "gray27", fill = "white") +
  xlab(expression(paste("carbon accumulation rate (t C ", ha^-1, yr^-1,")"))) +
  geom_point(aes(x = cmem_rates_cohort$acc_C[1], y = 0), color = colors[1], size = 3) +
  geom_point(aes(x = cmem_rates_cohort$acc_C[2], y = 0), color = colors[4], size = 3) +
  geom_point(aes(x = cmem_rates_cohort$acc_C[3], y = 0), color = colors[1], size = 3, shape = 17)+
  geom_point(aes(x = cmem_rates_cohort$acc_C[4], y = 0), color = colors[4], size = 3, shape = 17)+
  geom_segment(aes(x = cmem_rates_cohort$acc_C[1], y = 0, xend = cmem_rates_cohort$acc_C[1], yend = Inf),
               color = colors[1], size = 1.5, linetype = "dotted") +
  geom_segment(aes(x = cmem_rates_cohort$acc_C[2], y = 0, xend = cmem_rates_cohort$acc_C[2], yend = Inf),
               color = colors[4], size = 1.5, linetype = "dotted") + 
  geom_segment(aes(x = cmem_rates_cohort$acc_C[3], y = 0, xend = cmem_rates_cohort$acc_C[3], yend = Inf),
               color = colors[1], size = 1.5, linetype = "dotted") +
  geom_segment(aes(x = cmem_rates_cohort$acc_C[4], y = 0, xend = cmem_rates_cohort$acc_C[4], yend = Inf),
               color = colors[4], size = 1.5, linetype = "dotted") +
  theme_bw() +
  ylab("count") -> Fig4_panelC


## Figure 4d - surface elevation predictions for AGB only simulation ####

# Collect elevation values for each age-location cohort at the last year of the
# simulation
cmem_agb %>% 
  filter(iteration %in% c("corn-ancestral", "sellman-ancestral",
                          "corn-modern", "sellman-modern") & year == "2100") %>% 
  pull(value) -> end_points_agb

cmem_agb %>% 
  mutate(year = as.numeric(year)) %>%
  mutate(color_code = case_when(iteration == "corn-ancestral" ~ 1,
                                iteration == "sellman-ancestral" ~ 2,
                                iteration == "corn-modern" ~ 3,
                                iteration == "sellman-modern" ~ 4,
                                T ~ 0),
         size_code = case_when(iteration %in% c("corn-ancestral", "sellman-ancestral",
                                                "corn-modern", "sellman-modern") ~ "big",
                               T ~ "small")) %>%
  ggplot(aes(x = year, y = value, group = iteration, color = factor(color_code), size = factor(size_code))) +
  geom_line(aes(alpha = size_code)) +
  ylab("marsh elevation (cm NAVD88)") +
  scale_color_manual(values = c("gray11", colors[1], colors[4], colors[1], colors[4])) +
  scale_size_manual(values = c(1.5, 0.8)) +
  scale_alpha_manual(values = c(0.8, 0.1)) +
  geom_point(aes(x = 2100, y = end_points_agb[1]), color = colors[1], size = 3) +
  geom_point(aes(x = 2100, y = end_points_agb[2]), color = colors[4], size = 3) +
  geom_point(aes(x = 2100, y = end_points_agb[3]), color = colors[1], size = 3, shape = 17) +
  geom_point(aes(x = 2100, y = end_points_agb[4]), color = colors[4], size = 3, shape = 17) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(22.5,35) -> Fig4_panelD


## Figure 4e - compare variation for agb + bgb and abg only models ####

cmem_full %>% 
  filter(year == "2100") %>% 
  pull(value) -> surface_elevation_2100_full

cmem_agb %>% 
  filter(year == "2100") %>% 
  pull(value) -> surface_elevation_2100_agb

# Calculate decrease in variance from full to ag only
(var(surface_elevation_2100_full) - var(surface_elevation_2100_agb)) / var(surface_elevation_2100_full)

tibble(y = c(surface_elevation_2100_full, surface_elevation_2100_agb),
       x = rep(c("ag + bg", "ag only"), each = length(surface_elevation_2100_agb))) %>% 
  mutate(x = factor(x, levels = c("ag only", "ag + bg"))) %>% 
  ggplot(aes(x = x, y = y, fill = x)) +
  geom_violin(draw_quantiles = c(0.025,0.5, 0.975), size = 0.5, alpha = 0.4) +
  scale_fill_manual(values = c("gray67", "gray11")) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("elevation in year 2100 (cm NAVD88)") +
  xlab("scenario")-> Fig4_panelE

## Bring all plots together ####
Fig4_panelsBC <- plot_grid(Fig4_panelB, Fig4_panelC, nrow = 2, labels = c("b", "c"))
Fig4_panelA_label <- plot_grid(Fig4_panelA, labels = "a")
Fig4_panelsDE <- plot_grid(Fig4_panelD, Fig4_panelE, nrow = 1,
                                    labels = c("d", "e"), rel_widths = c(3,2))
Fig4_panelsABC <- plot_grid(Fig4_panelA_label, Fig4_panelsBC, rel_widths = c(3,2))
Fig4_panelsABCDE <- plot_grid(Fig4_panelsABC, Fig4_panelsDE, nrow = 2)

png(here("figs_tables", "Figure4.png"), height = 6.8, width = 8, units = "in", res = 300)
Fig4_panelsABCDE
dev.off()
