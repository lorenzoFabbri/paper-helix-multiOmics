# Script containing functions to perform analyses on HPC server
# Author: Lorenzo Fabbri

run_analysis <- function(choice, key.save.results) {
  source("code/ggm/ggm.R")
  
  save.data <- FALSE
  is.hpc <- FALSE
  
  ##### Parameters GGMs #####
  if (choice == "ggm") {
    path.save.res <- file.path(paste0("results/ggm/", key.save.results, "/"))
  } else {
    # Bootstrapping
    path.save.res <- file.path(paste0("results/ggm/", key.save.results, 
                                      "_boot/"))
  }
  active <- "corpcor"
  num.iters.boot <- 100
  params <- list(
    method.ggm = "prob", 
    cutoff.ggm = 0.8, 
    method.mixed = tidygraph::group_infomap, 
    path.save.plots = path.save.res
  )
  ##############################################################################
  
  if (choice == "ggm") {
    if (!dir.exists(path.save.res)) {
      dir.create(path = path.save.res)
    }
    
    ##### Pipeline #####
    ggms <- main.pipeline.ggm(params = list(
      omic.type = c("metabun", "metab_blood", "proteome", "methylome"), 
      package.corr = active, is.hpc = is.hpc, 
      key.save = key.save.results, 
      boot = list(
        perform = FALSE
      )
    ), save.data = save.data)
    if (save.data) {base::save(ggms, file = "../data/intermediate_res_ggm/ggms.RData")}
    processed.ggms <- process.ggms(list.ggms = ggms, active = active, 
                                   filter.mixed.interactions = FALSE, 
                                   is.directed = FALSE, 
                                   params = params, 
                                   boot = FALSE, save.data = save.data)
    if (save.data) {
      base::save(processed.ggms, 
                 file = "../data/intermediate_res_ggm/processed_ggms.RData")
    }
    is.stratified <- ggms[[names(ggms)[1]]]$params$stratification
    rm(ggms)
    gc()
    net.props <- net.properties(ggms = processed.ggms, type.networks = "ggm", 
                                path.save = path.save.res, is.merged = FALSE)
    if (is.stratified) {
      for (idx in 1:(length(processed.ggms) / 2)) {
        res <- merge.networks(ggms = processed.ggms, exposure.group = idx, 
                              how.to.join = how.to.join, save.data = save.data)
      }
    } else {
      res <- merge.networks(ggms = processed.ggms, exposure.group = 5, 
                            how.to.join = "inner", 
                            type.networks = "ggm", 
                            path.save = path.save.res, 
                            save.data = save.data)
      # Save name biomarkers for bootstrapping
      name.biomarkers <- res$mod_2.2.2.5.5$net %>%
        tidygraph::activate(what = "nodes") %>%
        tibble::as_tibble() %>%
        dplyr::mutate(name.ready = dplyr::case_when(
          label == "exposure" ~ paste0("log.", tolower(name), "_e"), 
          label == "proteome" ~ paste0(name, "_p"), 
          label == "serum metabolome" ~ paste0(name, "_ms"), 
          label == "urinary metabolome" ~ paste0(name, "_mu"), 
          label == "methylome" ~ paste0(name, "_me")
        ))
      if (save.data) {
        name.biomarkers %>%
        write.csv(file = paste0("data/merged_biomarkers_ggm.csv"), 
                  quote = FALSE, row.names = FALSE)
      }
    }
    if (save.data) {
      base::save(res, 
                 file = "../data/intermediate_res_ggm/merged_net.RData")
    }
    rm(processed.ggms)
    gc()
    net.merged.props <- net.properties(ggms = res, type.networks = "ggm", 
                                       path.save = paste0(path.save.res, "merged_"), 
                                       is.merged = TRUE, how.to.join = "inner")
    gc()
    
    analyse.cc(res$mod_2.2.2.5.5$net, path.save = path.save.res)
    # End if GGMs
    ############################################################################
    
  } else if(choice == "ggm.boot") { # Bootstrapping
    if (!dir.exists(path.save.res)) {
      dir.create(path = path.save.res)
    }
    
    for (itr in 1:num.iters.boot) {
      cat("\n####################################################")
      cat(paste0("\nBootstrapping. Iteration: ", itr, "...\n"))
      cat("####################################################\n")
      
      ggms <- main.pipeline.ggm(params = list(
        omic.type = c("metabun", "metab_blood", "proteome", "methylome"), 
        package.corr = active, is.hpc = is.hpc, 
        key.save = key.save.results, 
        boot = list(
          perform = TRUE, 
          seed = itr, 
          path.save.subjects = path.save.res
        )
      ), save.data = FALSE)
      processed.ggms <- process.ggms(list.ggms = ggms, active = active, 
                                     filter.mixed.interactions = FALSE, 
                                     is.directed = FALSE, 
                                     params = params, 
                                     boot = TRUE, save.data = save.data)
      rm(ggms)
      gc()
      res <- merge.networks(ggms = processed.ggms, exposure.group = 5, 
                            how.to.join = "inner", 
                            type.networks = "ggm", 
                            path.save = path.save.res, 
                            boot = TRUE, save.data = save.data)
      base::save(res, 
                 file = paste0(path.save.res, "merged_net", itr, ".RData"))
      rm(processed.ggms)
      gc()
    } # End loop iterations
    # End if bootstrapping
    ##########################################################################
      
    } else {
      stop(call. = TRUE)
    }
  
} # End function run analyses HPC

setwd("/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/")
options(device = NULL)
run_analysis(choice = "ggm", key.save.results = "tests")
