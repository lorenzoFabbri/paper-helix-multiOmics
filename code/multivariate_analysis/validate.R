# Author: Lorenzo Fabbri
# Script containing functions to perform validation by permutation

library(tibble)

source("./code/multivariate_analysis/model.R")

# Helper function to extract results from each run of permutation
extract.results <- function(res, iter) {
  
  measures <- res$perf$measures
  metrics  <- c("Q2", "Q2.total")
  
  # Iterate over relevant metrics
  for (metric.name in metrics) {
    metric  <- measures[[metric.name]]
    values  <- metric$values
    summary <- metric$summary
  }
}

# Main function to perform validation by permutation
validate.by.permutation <- function(data, pls.params, 
                                    num.iters.permut) {
  
  # Main loop to obtain `num.iters.permut` samples of each metric of interest
  res.list <- list()
  for (iter in 1:num.iters.permut) {
    cat("\n============================")
    cat(paste("\nIteration number: ", iter, "\n"))
    cat("============================\n")
    
    res.pls <- pls.simple(data, pls.params, validate = TRUE)
    # Append all results to list to process
    res.list[[paste0("iter_", iter)]] <- list(res.pls)
  }
  
  # Process list of results from all models
  for (iter in names(res.list)) {
    res <- res.list[[iter]]
    extracted.info <- extract.results(res, iter = iter)
  }
}
