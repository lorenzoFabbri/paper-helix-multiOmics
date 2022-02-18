# Author: Lorenzo Fabbri
# Script containing functions to pre-process and transform datasets

library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(labelled)
library(VIM)
library(imputeMissings)
library(tidyselect)
library(parallel)

source("code/multivariate_analysis/plot.R")

# Helper function to handle missing values
handle.missing.values <- function(data, method) {
  if (method == 'remove') {
    data <- data %>%
      tidyr::drop_na()
  } else if (method %in% c("kNN", "knn")) {
    data <- data %>%
      VIM::kNN(data = .) %>%
      dplyr::select(-tidyselect::ends_with("_imp"))
  } else if (method %in% c("median/mode")) {
    data <- data %>%
      imputeMissings::impute(., method = method)
  }
  
  return(data)
}

# Helper function to change type of given column
change.type <- function(data, col, new.type) {
  data[[col]] %<>% new.type
  
  return(data)
}

# Helper function to scale tibbles
scale.tibble <- function(data, ...) {
  
  if (... == "autoscale") {
    data <- (data - mean(data)) / sd(data)
    #data <- data / sd(data)
  } else if (... == "range") {
    data <- (data - mean(data)) / (max(data) - min(data))
    #data <- data / (max(data) - min(data))
  } else if (... == "identity") {
    data <- data
  } else {
    cat("\nMethod not yet implemented.\n")
    stop(call. = TRUE)
  }
  
  return(data)
}

# Helper function to scale datasets by column
scale.by.group <- function(data, group, start, end, method) {
  computing.version <- "parallel"
  
  if (group == "none") { # We do not scale data differently for each group
    
    if (computing.version == "parallel") {
      ##### Parallel version #####
      HelixID <- data$HelixID
      data.scaled <- parallel::mclapply(X = data %>%
                                          dplyr::select(-c(HelixID)), 
                                        FUN = scale.tibble, 
                                        ... = method, 
                                        mc.cores = parallel::detectCores() - 1) %>%
        tibble::as_tibble() %>%
        tibble::add_column(HelixID)
      #####                  #####
    } else {
      ##### Serial version #####
      data.scaled <- data %>%
        dplyr::mutate_at(vars(-c(HelixID)), 
                         .funs = list(scale.tibble), ... = method)
      #####                #####
    }
    
    # Add `group` anyway to plot by `group`
    group <- stringr::str_sub(data$HelixID, start = start, end = end)
    data.scaled <- tibble::add_column(data.scaled, group)
  } else { # We group by `group` and apply the scaling method to each sub-dataset
    group <- stringr::str_sub(data$HelixID, start = start, end = end)
    data <- tibble::add_column(data, group)
    
    data.scaled <- data %>%
      dplyr::group_by(group) %>%
      dplyr::mutate_at(vars(-c(HelixID, group)), 
                       .funs = list(scale.tibble), ... = method)
  }
  
  return(data.scaled)
}

# Main function to pre-process datasets
preprocess.data <- function(data, params) {
  validate <- params$validate
  
  # Drop rows containing missing values
  method.na <- params$method.na
  # cat("\nRemoving missing values...\n")
  # print(paste0("Dimension before removal: ", dim(data$exposures)[1], "."))
  # print(paste0("Dimension before removal: ", dim(data$omics)[1], "."))
  # data$exposures <- handle.missing.values(data$exposures, method = method.na)
  # data$omics     <- handle.missing.values(data$omics, method = method.na)
  # print(paste0("Dimension after removal: ", dim(data$exposures)[1], "."))
  # print(paste0("Dimension after removal: ", dim(data$omics)[1], "."))
  
  ret <- list(exposures = data$exposures, omics = data$omics, metadata = data$metadata)
  return(ret)
}

# Main function to transform (e.g., centering, scaling) datasets
transform.data <- function(data, scale.new, is.adjusted) {
  
  # Scale and center datasets by group
  if (scale.new$perform == TRUE) {
    cat("\nChecking if missing values...\n")
    if (sum(is.na(data$exposures)) != 0) {
      stop("Missing values in Exposome.")
    }
    if (sum(is.na(data$omics)) != 0) {
      stop("Missing values in Omics")
    }
    
    # For the Helix datasets, I assume the presence of a column
    # named `cohort` which contains the cohort ID
    # in a 3-letter format
    cat("\nScaling exposures...\n")
    data$exposures <- scale.by.group(data = data$exposures, group = scale.new$group, 
                                     start = 1, end = 3, method = scale.new$method)
    cat("\nScaling -omics...\n")
    # methylome.tmp <- data$omics %>%
    #   dplyr::select(dplyr::ends_with("_me"))
    data$omics <- scale.by.group(data = data$omics, #%>%
                                   #dplyr::select(-dplyr::ends_with("_me")), 
                                 group = scale.new$group, 
                                 start = 1, end = 3, method = scale.new$method)
    
    data$exposures <- data$exposures %>%
      dplyr::ungroup()
    data$omics <- data$omics %>%
      dplyr::ungroup()
    
    #data$omics <- dplyr::bind_cols(data$omics, methylome.tmp)
    #rm(methylome.tmp)
    gc()
  }
  
  # Drop columns containing metadata
  data$omics <- data$omics %>%
    dplyr::select(-c(HelixID))
  data$exposures <- data$exposures %>%
    dplyr::select(-c(HelixID))
  
  if (scale.new$perform == TRUE) {
    #boxplot.vars.by.group(dat = data$exposures, group = scale.new$group)
    #boxplot.vars.by.group(dat = data$omics, group = scale.new$group)
    
    data$exposures <- data$exposures %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(group))
    data$omics <- data$omics %>%
      dplyr::ungroup() %>%
      dplyr::select(-c(group))
  }
  
  ret <- list(exposures = data$exposures, omics = data$omics, metadata = data$metadata)
  return(ret)
}
