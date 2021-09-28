# ------------------------------- #
# Calculations of initial weights # 
# of propagules in pots           #
# ------------------------------- #

library(here)

# Read in data
pot_level <- read.csv(here("supp_data/", "InitialConditions.csv"))
weights <- read.csv(here("supp_data/", "PropaguleWeights.csv"))

head(pot_level)
head(weights)

# Create Id variable in weights dataset to match up
weights$Id <- as.character(interaction(weights$Site, weights$Layer, sep = ""))

# Create stem id variable for both datasets
weights$stemid <- as.character(interaction(weights$Id, weights$New.ID, sep = "_"))
pot_level$stemid <- as.character(interaction(pot_level$Id, pot_level$Number, sep = "_"))

# Number of unique stem ids: should be 384 (384 stems / 96 pots)
length(unique(pot_level$stemid))

merge(weights, pot_level, by = "stemid") %>% 
  select(pot = Pot,
         weight = Weight.x,
         genotype = Id.x) %>% 
  arrange(pot) -> ic_data

# Remove all other dataframes except the ones we need
rm(list= ls()[!(ls() %in% c('agb_adjust','bgb_adjusted_full', 'ic_data'))])
