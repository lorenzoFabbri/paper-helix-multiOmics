# Main script to run analyses and analyze results for reviews Env Int

rm(list=ls())
source("code/paper/images.R")
source("code/paper/tables.R")

################################################################################
##### Review 1
################################################################################

# Reviewer 4
## Supplementary tables for population description
tables.population.desc(to.save = TRUE, tbl.by = "Period")
tables.population.desc(to.save = TRUE, tbl.by = "Cohort")

