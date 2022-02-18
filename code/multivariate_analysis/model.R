# Author: Lorenzo Fabbri
# Script containing functions to fit PLS models

library(readr)
library(mixOmics)
library(parallel)
library(stringr)
library(dplyr)
library(fastDummies)
library(robust)
library(tibble)

source("code/multivariate_analysis/preprocess.R")
source("code/multivariate_analysis/plot.R")

# Function to perform 2-stage residual-outcome regression
perform.2sr <- function(outcome, covariates, method = "robust") {
    mc.cores <- parallel::detectCores() - 1
    
    ret1 <- parallel::mclapply(X = outcome %>%
                                   dplyr::select(tidyselect::ends_with(c("_p", "_mu"))), 
                               FUN = perform.2sr.single, 
                               covariates = covariates, is.serum = FALSE, 
                               method = method, 
                               mc.cores = mc.cores) %>%
        tibble::as_tibble()
    ret2 <- parallel::mclapply(X = outcome %>%
                                   dplyr::select(tidyselect::ends_with("_ms")), 
                               FUN = perform.2sr.single, 
                               covariates = covariates, is.serum = TRUE, 
                               method = method, 
                               mc.cores = mc.cores) %>%
        tibble::as_tibble()
    
    if (any(grepl("_me", colnames(outcome)))) {
        ret3 <- parallel::mclapply(X = outcome %>%
                                       dplyr::select(tidyselect::ends_with("_me")), 
                                   FUN = perform.2sr.single, 
                                   covariates = covariates, 
                                   is.serum = FALSE, is.methyl = TRUE, 
                                   method = method, 
                                   mc.cores = mc.cores) %>%
            tibble::as_tibble()
        ret <- dplyr::bind_cols(ret1, ret2, ret3)
        rm(list = c("ret1", "ret2", "ret3"))
    } else {
        ret <- dplyr::bind_cols(ret1, ret2)
        rm(list = c("ret1", "ret2"))
    }
    
    gc()

    return(ret)
}

######################################################
##### Helper function to perform 2sr on 1 column #####
######################################################
perform.2sr.single <- function(out, covariates, is.serum, is.methyl = FALSE, 
                               method) {
    # Here `out` is just one variable (i.e., one column)
    out <- tibble::as_tibble(out)
    covariates <- tibble::as_tibble(covariates)
    
    # Check if to be adjusted also for omic-specific confounders
    #if (grepl("_ms", names(out)) == FALSE) {
    if (is.serum == FALSE) {
        covariates <- covariates %>%
            dplyr::select(-c("hs_dift_mealblood_imp"))
    }
    if (is.methyl == TRUE) {
        # Methylation data we are using were adjusted for cohort, sex, age
        covariates <- covariates
            #dplyr::select(-c("e3_sex.x", "cohort.x", "age_sample_years.x"))
    }
    
    if (method == "standard") {
        res.model <- lm(formula = paste0(names(out), "~", 
                                         paste0(names(covariates), collapse = "+")), 
                        data = cbind(out, covariates))
    } else if (method == "robust") {
        res.model <- robust::lmRob(formula = paste0(names(out), "~", 
                                                    paste0(names(covariates), 
                                                           collapse = "+")), 
                                   data = cbind(out, covariates))
    }
    
    ret <- res.model[["residuals"]] + res.model[["coefficients"]][1]

    return(ret)
}

# Helper function to adjust for covariates
perform.adj <- function(data, metadata, covariates, method) {
    # Change variables type (chr to factor)
    metadata <- metadata %>%
        dplyr::mutate_at(dplyr::all_of(covariates$covariates.char), 
                         list(as.factor))
    
    covariates <- c(covariates$covariates.char, covariates$covariates.num)
    
    if (method == "residual") {
        data <- perform.2sr(outcome = data, covariates = metadata[, covariates])
    }
    else if (method == "variables") {
        # Simply concatenate the covariates to the matrix of descriptors (e.g., exposures)
        data <- data
    } else {
        cat("\nMethod not yet implemented.\n")
        stop(call. = TRUE)
    }

    return(data)
}

# Helper function to fit `mixOmics` PLS model
fit.model <- function(data, pls.params, no.sparsity, keepX, keepY) {
    if (pls.params$is.sparse == FALSE) {
        model <- mixOmics::pls(X = data$exposures, Y = data$omics, 
                               ncomp = pls.params$num.components, 
                               scale = FALSE, mode = "regression", 
                               max.iter = pls.params$max.iter, near.zero.var = FALSE)
    } else {
        if (no.sparsity == TRUE) {
            model <- mixOmics::spls(X = data$exposures, Y = data$omics, ncomp = pls.params$num.components,
                mode = "regression", scale = FALSE, max.iter = pls.params$max.iter)
        } else {
            model <- mixOmics::spls(X = data$exposures, Y = data$omics, ncomp = pls.params$num.components,
                mode = "regression", scale = FALSE, max.iter = pls.params$max.iter,
                keepX = keepX, keepY = keepY)
        }
    }

    return(model)
}

# Helper function to run main pipeline (fit, tune, perf from `mixOmics`)
run.pipeline <- function(data, pls.params, model.code) {

    # Fit initial model
    model <- fit.model(data = data, pls.params = pls.params, no.sparsity = TRUE,
        keepX = NULL, keepY = NULL)

    # Check performances
    perf.obj <- mixOmics::perf(object = model, validation = "Mfold", folds = pls.params$nfolds, 
                               nrepeat = pls.params$nrepeat, criterion = "all")

    # Tune hyper-parameters
    if (pls.params$is.sparse == TRUE) {
        # If we adjust the models for some confounders, we must keep all the 
        # variables in the X matrix
        pls.params$keepX <- c(3, 5, 7, 10, ncol(data$exposures))
        pls.params$keepY <- c(3, 5, 7, 10, ncol(data$omics))
        
        tune.obj <- mixOmics::tune(method = "spls", 
                                   X = data$exposures, Y = data$omics, 
                                   ncomp = pls.params$num.components, 
                                   test.keepX = pls.params$keepX, 
                                   test.keepY = pls.params$keepY, 
                                   mode = "regression", 
                                   nrepeat = pls.params$nrepeat, validation = "Mfold", folds = pls.params$nfolds, 
                                   measure = "cor", 
                                   center = FALSE, scale = FALSE, 
                                   max.iter = pls.params$max.iter, 
                                   cpus = parallel::detectCores())
        print(plot(tune.obj, measure = "cor"))

        # Check again performances with re-fitted model
        #pls.params$num.components <- tune.obj$choice.ncomp
        
        ##### ATTENTION! #####
        # WE FORCE THE MODEL TO USE 3 VARIATES (COMPONENTS) FOR BOTH PLS AND 
        # sPLS, NO MATTER THE RESULTS OF `tune`!
        ######################
        
        model <- fit.model(data = data, 
                           pls.params = pls.params, no.sparsity = FALSE, 
                           keepX = tune.obj$choice.keepX, 
                           keepY = tune.obj$choice.keepY)
        perf.obj <- mixOmics::perf(object = model, validation = "Mfold", folds = pls.params$nfolds, 
                                   nrepeat = pls.params$nrepeat, criterion = "all")
        print(plot(perf.obj, criterion = "Q2.total"))
        
        cat("\n==================================================\n")
        cat("Dimension of Exposome (X):\n")
        print(dim(data$exposures))
        cat("Dimension of Omics (Y):\n")
        print(dim(data$omics))
        cat("\n--------------------------------------------------\n")
        cat("Summary of (s)PLS\n")
        cat("Number of components selected:\n")
        print(tune.obj$choice.ncomp)
        cat("Number of keepX selected:\n")
        print(tune.obj$choice.keepX)
        cat("Number of keepY selected:\n")
        print(tune.obj$choice.keepY)
        cat("\n==================================================\n")
        
        # Write information to file
        write(dim(data$exposures), 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        write(dim(data$omics), 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        write(tune.obj$choice.ncomp, 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        write(tune.obj$choice.keepX, 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        write(tune.obj$choice.keepY, 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        
    } else {
        tune.obj <- NULL
        cat("\n==================================================\n")
        cat("Dimension of Exposome (X):\n")
        print(dim(data$exposures))
        cat("Dimension of Omics (Y):\n")
        print(dim(data$omics))
        cat("\n==================================================\n")
        
        # Write information to file
        write(dim(data$exposures), 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
        write(dim(data$omics), 
              file = paste0(pls.params$path.save, 
                            "tables/", model.code, ".txt", collapse = ""), 
              append = TRUE)
    }
    
    # Plot predicted vs. actual values of -omics for all components
    plot.predicted(model)

    ret <- list(model = model, tune = tune.obj, perf = perf.obj, data = data)
    return(ret)
}

# Main function to fit PLS models
pls.simple <- function(data, pls.params, validate, covariates) {
    model.code <- pls.params$model.code
    
    # Pre-process datasets
    data <- preprocess.data(data = data, validate = validate)

    # Adjust for covariates using 2-stage regression
    if (pls.params$perform.adj$perform == TRUE) {
        
        if (pls.params$perform.adj$method == "variables") {
            stop(call. = TRUE)
            
            # Transform to factors, create dummy variables and append to 
            # matrix of descriptive variables (e.g., exposures)
            old.nrow <- nrow(ret$exposures)
            
            # ret$metadata <- ret$metadata %>%
            #     dplyr::mutate_at(dplyr::all_of(covariates$covariates.char), 
            #                      list(as.factor))
            metadata.tmp <- fastDummies::dummy_cols(.data = ret$metadata[, covariates$covariates.char], 
                                                    select_columns = covariates$covariates.char, 
                                                    remove_first_dummy = FALSE, 
                                                    remove_selected_columns = TRUE)
            ret$metadata <- cbind(metadata.tmp, ret$metadata[, covariates$covariates.num])
            # Add labels to confounders
            for (col in colnames(ret$metadata)) {
                labelled::var_label(ret$metadata[[col]]) <- "confounders"
            }
            
            # Combine exposures (actual predictors) w/ confounders
            ret$exposures <- cbind(ret$exposures, ret$metadata)
            
            if (old.nrow != nrow(ret$exposures)) {stop(call. = TRUE)}
        } else {
            omics <- perform.adj(data = data$omics %>% dplyr::select(-c(SampleID, HelixID)), 
                                 metadata = data$metadata, 
                                 covariates = covariates, 
                                 method = pls.params$perform.adj$method)
            data$omics <- base::cbind(data$omics[, c("SampleID", "HelixID")], 
                                      omics)
        }
    }
    
    # Transform datasets
    data <- transform.data(data = data, scale.new = pls.params$scale.new, 
                           is.adjusted = pls.params$perform.adj$perform)
    
    # Fit model
    ret <- run.pipeline(data = data, pls.params = pls.params, 
                        model.code = model.code)
    
    return(ret)
}
