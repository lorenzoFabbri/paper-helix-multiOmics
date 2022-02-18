# Script to produce preliminary results with sgPLS and Multi-Omics data
# Author: Lorenzo Fabbri

rm(list = ls())
library(sgPLS)

source("code/comparison_methods.R")
load("code/env")

##############################
##### DATA PREPROCESSING #####
##############################

covariates.char <- c("e3_sex.x")
covariates.num  <- c("age_sample_years.x", "zBMI")
covariates <- list(covariates.char = covariates.char, 
                   covariates.num = covariates.num)
omic.type <- "metabun"

params = list(
  time.point = 1, 
  omic.type = omic.type, 
  scale.new = list(perform = TRUE, 
                   group = "cohort", method = "range"), 
  validate = FALSE, 
  adjust = list(perform = TRUE, method = "residual", list = covariates)
)

dat <- generate.clean.data(params = params)

# Grouping for sgPLS
exposures <- list(phenols = c("mepa", "etpa", "prpa", "bupa", 
                              "bpa", "trcs", "oxbe"), 
                  phthalates = c("mep", "mibp", "mnbp", "mbzp", "mehp", 
                                "mehhp", "meohp", "mecpp", "ohminp", "oxominp"), 
                  pesticides = c("dmtp", "dep", "detp"))
print(length(unlist(exposures)))

# Short exposure names
dat$data.preprocessed$exposures <- dat$data.preprocessed$exposures %>%
  dplyr::rename_with(., ~ gsub("log.hcp_", "", .x)) %>%
  dplyr::rename_with(., ~ gsub("_cadj", "", .x))

# Order columns exposures according to list
dat$data.preprocessed$exposures <- dat$data.preprocessed$exposures[unlist(exposures)]
exposure.groups <- c(length(exposures$phenols), 
                     length(exposures$phenols) + length(exposures$phtalates), 
                     length(exposures$phenols) + length(exposures$phtalates) + length(exposures$pesticeds))

####################
##### MODELING #####
####################

X <- dat$data.preprocessed$exposures
Y <- dat$data.preprocessed$omics

ncomp <- 5
grid.X <- c(1, 3, ncol(X) / 2, 15, ncol(X))

##### sPLS
# Iterative strategy
spls.tune <- sgPLS::tuning.sPLS.X(X = X, Y = Y, 
                                  folds = 5, validation = "Mfold", 
                                  ncomp = 1, 
                                  keepX = NULL, grid.X = grid.X, 
                                  setseed = 1)

##### gPLS
gpls.fit <- sgPLS::gPLS(X = X, Y = Y, 
                        ncomp = ncomp, mode = "regression", 
                        keepX = rep(3, ncomp), #keepY = , 
                        ind.block.x = exposure.groups, #ind.block.y = , 
                        scale = FALSE, 
                        max.iter = 500, tol = 1e-06)

# Iterative strategy
gpls.tune <- sgPLS::tuning.gPLS.X(X = X, Y = Y, 
                                  folds = 5, validation = "Mfolds", 
                                  ncomp = 1, keepX = NULL, 
                                  grid.X = grid.X, ind.block.x = exposure.groups, 
                                  setseed = 1)

##### sgPLS
