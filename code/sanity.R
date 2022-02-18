# Author: Lorenzo Fabbri
# Script containing functions to perform sanity checks on obtained results 
# (e.g., compare results of custom script and package function)

library(mixOmics)

####################
##### mixOmics #####
####################

##### Helper function to load results mixOmics and data
load.res.mixomics <- function(model.code) {
  
  path.models <- "./results/models/images/"
  model.code <- paste0("mod_", model.code)
  path.model <- paste0(path.models, model.code)
  
  # Load parameters PLS
  params <- get(load(paste0(path.model, "/pls.params.Rdata")))
  
  # Load results mixOmics
  res <- get(load(paste0(path.model, "/ret.Rdata")))
  
  # Load data
  dat <- get(load(paste0(path.model, "/data.preprocessed.Rdata")))
  
  return(list(
    fitted.model = res, 
    data = dat, 
    params = params
  ))
  
}

##### Function to plot results PLS from mixOmics #####
check.pls.mixomics <- function(model.code, components) {
  
  loaded <- load.res.mixomics(model.code)
  res <- loaded$res
  dat <- loaded$data
  
  # Plot results using mixOmics's functions
  mixOmics::plotVar(object = res$model, comp = components)
  mixOmics::plotLoadings(object = res$model, block = "X", comp = components[1])
  mixOmics::plotLoadings(object = res$model, block = "X", comp = components[2])
  mixOmics::plotLoadings(object = res$model, block = "Y", comp = components[1])
  mixOmics::plotLoadings(object = res$model, block = "Y", comp = components[2])
  
  mixOmics::network(mat = res$model, comp = components)
  
  # Tables
  print(res$model$prop_expl_var)
  
}

##### Function to compare Relevance Networks with those of xMWAS
check.rn.xmwas <- function(model.code, num.components) {
  
  loaded <- load.res.mixomics(model.code)
  params <- loaded$params
  res <- loaded$res
  dat <- loaded$data
  
}

################
##### GGMs #####
################
