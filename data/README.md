The metadata for the XX data files in this folder are as follows:

## BelowgroundEvolution_PotLevel.csv
* **pot** = pot number. Numbers correspond to overall treatments/groupings described in detail below: 1-18 = ancestral polycultures, 19-30 = modern polycultures, 31-48 = mixed polycultures, 49-72 = modern monocultures, 73-96 = ancestral monocultures. 
* **frame** = which of the four frames the experimental plot (PVC) was located. Ranges from 1-4, with 1 being the closest to the main boardwalk at GCREW and 4 being the furthest.
* **depth** = distance between top of pot and soil level at time of harvest (cm). 
* **ln_depth** = natural log of the distance between top of pot and soil level at time of harvest (cm). 
* **ic_weight** = initial summed wet weight of the four propagules within each pot (g).
* **diversity** = diversity treatment ('mono' = monoculture, 'poly' = polyculture). Monocultures mean that four propagules of the same genotype were used to establish the experimental plot. Polycultures mean that four propagules of four different genotypes were used to establish the experimental plot.
* **age** = age cohort (ancestral, modern, mix). Genotypes were grouped into two age cohorts based on the depth that the seed that resulted in the germinated plant was collected from. Here, ancestral genotypes are dated between 1940 and 1980 and modern genotypes are dated between 2000 and 2020. A mix of two modern and two ancestral genotypes was used for some replicates ('mixed'). 
* **provenance** = site location where core for seed extraction was collected (Sellman, Hog, Kirkpatrick, Corn; all within the property of the Smithsonian Environmental Research Center in Edgewater, MD).
* **genotype** = unique genotype ID for each of the 16 genotypes. The first letter of each genotype reflects the provenance (*e.g.* C = Corn) and the last letter reflects the depth within the core, where A is the deepest layer of the collected core.
* **agb** = dried aboveground biomass (g).
* **bgb** = dried belowground biomass (g; includes roots and rhizomes).
* **rs** = root-to-shoot ratio (g/g; technically the ratio of roots *and* rhizomes to shoots; **bgb**/**agb**).
* **mean_tot_height** = mean height of stems (cm) at the time of harvest.
* **mean_mid_width** = mean width of stems (mm) at the time of harvest. Widths were taken within the middle third of the total stem height.
* **density** = live stem density.
* **beta** = parameter describing the rooting depth distribution (see supp_code/RootDepth.R for details). Higher values indicate more belowground biomass is allocated to greater depths within the soil.
