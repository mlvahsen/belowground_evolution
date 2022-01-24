# Load libraries
library(tidyverse); library(cowplot); library(here); library(geomtextpath)

cmem_full <- read_rds(here("outputs/CMEM_runs", "CMEM_predictions_full.rds"))

cmem_cohort <- cmem_full %>% filter(iteration %in% c("corn-ancestral", "corn-modern")) %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(cohort = case_when(iteration == "corn-ancestral" ~ "ancestral cohort",
                            T ~ "descendant cohort"))

cmem_full %>% 
  mutate(year = as.numeric(year)) %>%
  filter(!is.na(as.numeric(iteration))) %>% 
  ggplot(aes(x = year, y = value, group = iteration)) +
  geom_line(color = "gray11", alpha = 0.1, size = 0.8) +
  geom_textline(data = cmem_cohort, aes(x = year, y = value, label = cohort),
                size = 4, hjust = 0.9, linewidth = 2, color = rep(c("orange", "dodgerblue"), each = 80),
                fontface = 2)+
  ylab("marsh elevation (cm NAVD88)") +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(22,36) -> Fig4_evoonly

ggsave(here("figs_tables", "fig4_evoonly.png"), height = 2.8, width = 4.2, units = "in")
