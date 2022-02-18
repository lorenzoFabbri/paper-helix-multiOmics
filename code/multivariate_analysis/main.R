# Author: Lorenzo Fabbri
# Main script to run Multivariate Analyses

source("./code/multivariate_analysis/model.R")
source("./code/multivariate_analysis/plot.R")
source("./code/multivariate_analysis/tables.R")
source("./code/multivariate_analysis/dictionaries.R")

library(rmarkdown)
library(purrr)
library(stringr)
library(labelled)

perform.analysis.main <- function(data, 
                                  type.pls, pls.params, validate) {
  
  # Covariates for adjusting models
  covariates <- list.covariates()
  
  # Stratification by exposure family/group
  exposures <- dict.exposure.groups()
  if (pls.params$stratification) {
    exposures.long <- list()
    for (grp in names(exposures)) {
      exposures.long[[grp]] <- sapply(exposures[[grp]], function(x) {
        paste0("log.hcp_", x, "_cadj")
      }) %>% unname()
    }
    
    data$exposures <- data$exposures %>%
      dplyr::select(., dplyr::one_of(c("HelixID", "SubjectID", 
                                       exposures.long[[pls.params$stratification.group]])))
  } # End stratification
  
  # Write information to file
  for (el in c(pls.params$model.code, 
               pls.params$omic.type, pls.params$time.point, 
               pls.params$is.sparse, 
               pls.params$perform.adj$method, 
               pls.params$scale.new$method, pls.params$scale.new$group, 
               pls.params$max.iter, pls.params$nrepeat, pls.params$nfolds)) {
    write(el, 
          file = paste0("./results/models/tables/", 
                        pls.params$model.code, ".txt", collapse = ""), 
          append = TRUE)
  }
  
  # Select and run analyses
  if (type.pls == "simple") {
    res.pls <- pls.simple(data = data, 
                          pls.params = pls.params, validate = validate, 
                          covariates = covariates)
    model <- res.pls$model
    tune  <- res.pls$tune
    perf  <- res.pls$perf
    data.preprocessed <- res.pls$data
  }
  
  if (validate == FALSE) {
    # Transfer labels exposures to processed data
    for (col in colnames(data.preprocessed$exposures)) {
      if (col %in% colnames(data$exposures)) {
        lab.tmp <- labelled::var_label(data$exposures[[col]])
      } else { # Not really safe but still
        lab.tmp <- "confounders"
      }
      
      labelled::var_label(data.preprocessed$exposures[[col]]) <- lab.tmp
    }
    
    # Plot results
    if (pls.params$plotting$perform == TRUE) {
      plot.internal(model = model, perf.obj = perf, 
                    metadata = data.preprocessed$metadata, 
                    pls.params = pls.params)
      plot.perf(model = model, perf.obj = perf, 
                data = data.preprocessed, 
                pls.params = pls.params)
    }
  }
  
  ret <- list(model = model, tune = tune, perf = perf)
  # Save selected objects to file
  obj.to.save <- c("data.preprocessed", "pls.params", "ret")
  dir.create(file.path("./results/models/images", pls.params$model.code), 
             showWarnings = FALSE)
  for (obj in obj.to.save) {
    save(list = obj, file = paste0("./results/models/images/", 
                                   pls.params$model.code, "/", obj, ".Rdata"))
  }
  
  return(ret)
  
}

# Function to render R Notebook with arguments
# render.report <- function(notebook.params, ...) {
#   template   <- "./code/multivariate_analysis.Rmd"
#   out.format <- "html_document"
#   
#   rmarkdown::render(template, 
#                     output_dir = "./results/models", 
#                     output_file = notebook.params$model.code, 
#                     output_format = out.format, 
#                     params = notebook.params)
#   invisible(TRUE)
#   
# }

render.report <- function(params) {
  
  # Load all data: -omics, exposures, metadata
  load("./code/env")
  
  # Parameters
  time.point <- ifelse(params$time.point == 1, "", "2")
  exposome <- get(paste0("exp_", params$omic.type, time.point))
  omic     <- get(paste0("omic_", params$omic.type, time.point))
  metadata     <- get(paste0("metadata", time.point))
  metadata_all <- readr::read_csv("./data/metadata_all.csv", col_names = TRUE)
  metadata     <- metadata_all %>%
    filter(period == ifelse(params$time.point == 1, "A", "B")) %>%
    inner_join(metadata, by = "HelixID")
  
  dat <- list(
    exposures = exposome, 
    omics     = omic, 
    metadata  = metadata
  )
  
  # Run main analysis
  ret <- perform.analysis.main(data = dat, 
                               type.pls = "simple", pls.params = params, 
                               validate = FALSE)
  
  return(ret)
}

# Function to loop over all required arguments
render.reports <- function() {
  rm(list = ls())
  
  # Prepare directories
  if (file.exists(file.path("./results/models"))) {
    unlink(file.path("./results/models"), recursive = TRUE)
    dir.create(file.path("./results/models"), showWarnings = FALSE)
  } else {
    dir.create(file.path("./results/models"), showWarnings = FALSE)
  }
  dir.create(file.path("./results/models", "images"), 
             showWarnings = FALSE)
  dir.create(file.path("./results/models", "plots_by_time"), 
             showWarnings = FALSE)
  dir.create(file.path("./results/models", "tables"), 
             showWarnings = FALSE)
  
  # Dictionary to create correct model.code
  codes <- params.to.model.code(model.type = "xPLS")
  code.scaling.type <- codes$code.scaling.type
  code.adj.type <- codes$code.adj.type
  code.model.type <- codes$code.model.type
  #code.time.type
  code.omic.type <- codes$code.omic.type
  code.exposure.type <- codes$code.exposure.type
  
  # Parameters
  omic.type <- names(codes$code.omic.type)
  stratification.group <- names(codes$code.exposure.type)
  stratification <- c(TRUE)
  time.point <- c(1, 2)
  num.components.pls <- 5
  num.components.spls <- 5
  is.sparse <- c(FALSE)
  perform.adj <- list(perform = TRUE, method = "residual")
  log.transform <- FALSE
  scale.new <- list(perform = TRUE, group = "none", 
                    method = "autoscale")
  max.iter <- 500
  nrepeat  <- 25
  nfolds   <- 5
  color.by.indiv <- "zBMI"
  plotting <- list(perform = TRUE, save = TRUE)
  path.save <- "./results/models/"
  
  # Create dataframe with all possible combinations of parameters
  all.params <- list(omic.type = omic.type, 
                     stratification.group = stratification.group, 
                     stratification = stratification, 
                     time.point = time.point, 
                     is.sparse = is.sparse, 
                     log.transform = log.transform, 
                     max.iter = max.iter, nrepeat = nrepeat, nfolds = nfolds, 
                     color.by.indiv = color.by.indiv)
  all.params <- expand.grid(all.params)
  print(all.params)
  
  # Iterate over combinations (rows) of parameters and perform analysis
  for (row in 1:nrow(all.params)) {
    el <- list()
    for (param in names(all.params)) {
      el[[param]] <- all.params[row, param]
    }
    
    el$scale.new   <- scale.new
    el$perform.adj <- perform.adj
    el$plotting    <- plotting
    
    el$num.components <- ifelse(el$is.sparse == TRUE, 
                                num.components.spls, num.components.pls)
    model <- ifelse(el$is.sparse == TRUE, "spls", "pls")
    el$model.code <- paste0("mod_", 
                            code.scaling.type[[el$scale.new$method %>% as.character()]], 
                            code.adj.type[[el$perform.adj$method %>% as.character()]], 
                            code.model.type[[model %>% as.character()]], 
                            el$time.point, 
                            code.omic.type[[el$omic.type %>% as.character()]], 
                            code.exposure.type[[el$stratification.group %>% as.character()]])
    el$path.save <- path.save
    
    # Create text file to save results analysis
    file.name <- paste0("./results/models/tables/", 
                        el$model.code, ".txt")
    file.create(file.name)
    
    # Create report
    render.report(el)
  }
  
  # Create tables results and grid of plots by time point
  if (el$plotting$perform & el$plotting$save) {
    make.report.by.time(path.all.plots = "./results/models/plots_by_time/", 
                        pcs = el$num.components)
  }
  
}
