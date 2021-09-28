# Remove stem weight from belowground layers 0-5cm and 5-10cm for GCREW
# diversity 2018 analysis

library(tidyverse); library(here)

bgb_og <- read_csv(here("supp_data", "BelowgroundBiomass.csv"))
bgb_stem <- read_csv(here("supp_data", "BelowgroundRemoveStems.csv"))

# Filter out data that is 0-5 for all pots
bgb_og %>% 
  filter(depth_roots %in% c(0,5)) %>% 
  mutate(depth_roots_lower = depth_roots + layer_width) %>% 
  arrange(pot, depth_roots, depth_roots_lower) -> bg_og_05 

# Fix pot 39 (the one that half was powerwashed) in stems data sheet to match up
# with original data sheet
bgb_stem %>% 
  mutate(depth_roots_lower = ifelse(depth_roots == "0-15", 15, as.numeric(depth_roots) + 5)) %>% 
  mutate(depth_roots = str_extract(depth_roots, "[0-9]+")) %>% 
  mutate(depth_roots = as.numeric(depth_roots)) %>% 
  add_row(pot = 39, depth_roots = 5, biomass_type = "stem", "biomass" = 0, notes = "check", depth_roots_lower = 10) %>% 
  arrange(pot, depth_roots, depth_roots_lower)-> bgb_stem
# Warning does not seem to mess up anything

# Check alignment
bgb_stem$pot - bg_og_05$pot
bgb_stem$depth_roots - bg_og_05$depth_roots
bgb_stem$depth_roots_lower - bg_og_05$depth_roots_lower
# Should all be zeros

# Create new weight columns with and without stem weights
bg_og_05 %>% 
  mutate(biomass_stems = bgb_stem$biomass,
         biomass_nostems = biomass - bgb_stem$biomass) -> bgb_adjusted

# Join in with the rest of the bgb data
bgb_og %>% 
  filter(depth_roots > 5) %>% 
  mutate(depth_roots_lower = depth_roots + layer_width,
         biomass_stems = 0,
         biomass_nostems = biomass) -> bgb_og_rest

rbind(bgb_adjusted, bgb_og_rest) %>% 
  arrange(pot, depth_roots, depth_roots_lower) -> bgb_adjusted_full

# Figure out extra stem weight that will need to be added to aboveground biomass
bgb_adjusted_full %>% 
  group_by(pot) %>% 
  summarize(stem_agb = sum(biomass_stems)) -> agb_adjust

# Remove all other dataframes except the ones we need
rm(list= ls()[!(ls() %in% c('agb_adjust','bgb_adjusted_full'))])
