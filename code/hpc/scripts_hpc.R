# Script containing functions to perform analyses on HPC server
# Author: Lorenzo Fabbri

run_analysis <- function(choice, key.save.results) {
  source("code/ggm/ggm.R")
  
  ##### Parameters GGMs #####
  if (choice == "ggm") {
    path.save.res <- file.path(paste0("results/ggm/", key.save.results, "/"))
  } else {
    # Bootstrapping
    path.save.res <- file.path(paste0("results/ggm/", key.save.results, 
                                      "_boot/"))
  }
  active <- "corpcor"
  num.iters.boot <- 50
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
      package.corr = active, is.hpc = TRUE, 
      key.save = key.save.results, 
      boot = list(
        perform = FALSE
      )
    ))
    base::save(ggms, file = "../data/intermediate_res_ggm/ggms.RData")
    cat("Processing models...\n")
    processed.ggms <- process.ggms(list.ggms = ggms, active = active, 
                                   filter.mixed.interactions = FALSE, 
                                   is.directed = FALSE, 
                                   params = params)
    base::save(processed.ggms, 
               file = "../data/intermediate_res_ggm/processed_ggms.RData")
    is.stratified <- ggms[[names(ggms)[1]]]$params$stratification
    rm(ggms)
    gc()
    cat("Computing network properties...\n")
    net.props <- net.properties(ggms = processed.ggms, type.networks = "ggm", 
                                path.save = path.save.res, is.merged = FALSE)
    if (is.stratified) {
      for (idx in 1:(length(processed.ggms) / 2)) {
        res <- merge.networks(ggms = processed.ggms, exposure.group = idx, 
                              how.to.join = how.to.join)
      }
    } else {
      cat("Merging networks...\n")
      res <- merge.networks(ggms = processed.ggms, exposure.group = 5, 
                            how.to.join = "inner", 
                            type.networks = "ggm", 
                            path.save = path.save.res)
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
        )) %>%
        write.csv(file = paste0("data/merged_biomarkers_ggm.csv"), 
                  quote = FALSE, row.names = FALSE)
    }
    base::save(res, 
               file = "../data/intermediate_res_ggm/merged_net.RData")
    cat("Computing network properties...\n")
    rm(processed.ggms)
    gc()
    net.merged.props <- net.properties(ggms = res, type.networks = "ggm", 
                                       path.save = paste0(path.save.res, "merged_"), 
                                       is.merged = TRUE, how.to.join = "inner")
    gc()
    
    cat("Analysing connected components...\n")
    analyse.cc(res$mod_2.2.2.5.5$net, path.save = path.save.res)
    # End if GGMs
    ############################################################################
    
  } else if(choice == "ggm.boot") { # Bootstrapping
    if (!dir.exists(path.save.res)) {
      dir.create(path = path.save.res)
    }
    
    for (itr in 1:num.iters.boot) {
      cat("\n####################################################")
      cat(paste0("\nBootstrapping. Iteration: ", itr, "...\n\n"))
      cat("####################################################\n")
      
      ggms <- main.pipeline.ggm(params = list(
        omic.type = c("metabun", "metab_blood", "proteome", "methylome"), 
        package.corr = active, is.hpc = TRUE, 
        key.save = key.save.results, 
        boot = list(
          perform = TRUE, 
          seed = itr
        )
      ))
      processed.ggms <- process.ggms(list.ggms = ggms, active = active, 
                                     filter.mixed.interactions = FALSE, 
                                     is.directed = FALSE, 
                                     params = params)
      rm(ggms)
      gc()
      res <- merge.networks(ggms = processed.ggms, exposure.group = 5, 
                            how.to.join = "inner", 
                            type.networks = "ggm", 
                            path.save = path.save.res)
      base::save(res, 
                 file = paste0(path.save.res, "merged_net", itr, ".RData"))
      rm(processed.ggms)
      gc()
    } # End loop iterations
    # End if bootstrapping
    ##########################################################################
      
    } else if (choice == "ggm.boot.par") {
      if (!dir.exists(path.save.res)) {
        dir.create(path = path.save.res)
      }
      
      ret <- bettermc::mclapply(1:num.iters.boot, function(itr) {
        
        ggms <- main.pipeline.ggm(params = list(
          omic.type = c("metabun", "metab_blood", "proteome", "methylome"), 
          package.corr = active, is.hpc = TRUE, 
          key.save = key.save.results, 
          boot = list(
            perform = TRUE, 
            seed = itr
          )
        ))
        processed.ggms <- process.ggms(list.ggms = ggms, active = active, 
                                       filter.mixed.interactions = FALSE, 
                                       is.directed = FALSE, 
                                       params = params)
        rm(ggms)
        gc()
        res <- merge.networks(ggms = processed.ggms, exposure.group = 5, 
                              how.to.join = "inner", 
                              type.networks = "ggm", 
                              path.save = path.save.res)
        base::save(res, 
                   file = paste0(path.save.res, "merged_net", itr, ".RData"))
        rm(processed.ggms)
        gc()
      }, mc.cores = 2)
      
      # End if parallel bootstrapping
      ##########################################################################
    }
  
} # End function run analyses HPC

setwd("/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/")
options(device = NULL)
run_analysis(choice = "ggm.boot", key.save.results = "debug")
