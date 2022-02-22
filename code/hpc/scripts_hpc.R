# Script containing functions to perform analyses on HPC server
# Author: Lorenzo Fabbri

run_analysis <- function(choice, key.save.results) {
  
  if (choice == "ggm") {
    source("code/ggm/ggm.R")
    
    ##### Parameters #####
    active <- "corpcor"
    path.save.res <- file.path(paste0("results/ggm/", key.save.results, "/"))
    if (!dir.exists(path.save.res)) {
      dir.create(path = path.save.res)
    }
    
    params <- list(
      method.ggm = "prob", 
      cutoff.ggm = 0.8, 
      method.mixed = tidygraph::group_infomap, 
      path.save.plots = path.save.res, 
      boot = list(
        perform = FALSE
      )
    )
    
    ##### Pipeline #####
    ggms <- main.pipeline.ggm(params = list(
      omic.type = c("metabun", "metab_blood", "proteome", "methylome"), 
      package.corr = active, is.hpc = FALSE, 
      key.save = key.save.results
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
        write.csv(file = paste0(path.save.res, "merged_biomarkers.csv"), 
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
  } # End if GGMs
  
} # End function run analyses HPC

setwd("/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/")
options(device = NULL)
run_analysis(choice = "ggm", key.save.results = "debug")
