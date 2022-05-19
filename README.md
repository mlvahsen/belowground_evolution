# Belowground Evolution Manuscript
Code and data for *Schoenoplectus americanus* belowground evolution manuscript

All code was run in R version 4.0.3 (2020-10-10). Necessary packages to run each script are identified at the top of each script under the heading "Load libaries". 

All R scripts are written assuming the below organizational structure:
```
belowground_evolution
│   README.md
└─── data
└─── figs_tables
└─── main_code
└─── outputs
└─── supp_code
└─── supp_data
```
The main script to run the hierarchical models is BelowgroundEvolution_script.R. This script needs to be run to generate the coda objects (*i.e.* generated samples from the Bayesian models) necessary to create figures 2-4 and S3-S4.  

Code to generate reported means and credible intervals in text is at the end of the code for the related figure. 

