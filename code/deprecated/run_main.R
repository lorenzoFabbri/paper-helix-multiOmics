# Temporary script to test new code for multivariate analysis

library(magrittr)

rm(list = ls())
load("./code/env")

change.type <- function(data, col, new.type) {
  data[[col]] %<>% new.type
  return(data)
}

exposome <- change.type(exp_metabun, "SubjectID", as.character)
omic     <- change.type(omic_metabun, "SubjectID", as.character)
metadata <- change.type(metadata, "SubjectID", as.character)
dat <- list(
  exposures = exposome, 
  omics = omic, 
  metadata = metadata
)

source("./code/multivariate_analysis/main.R")
pls.params <- list(
  num.components = 7, 
  num.vars       = 5, 
  is.sparse      = TRUE, 
  include.bmi    = FALSE, 
  perform.adj    = FALSE, 
  log.transform  = FALSE, 
  scale.new      = list(perform = TRUE, group = "cohort"), 
  max.iter       = 100, 
  nrepeat = 1, nfolds = 5
)
ret <- perform.analysis.main(dat, "simple", pls.params)
