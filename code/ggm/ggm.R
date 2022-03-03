# Script containing functions to generate GGMs
# Author: Lorenzo Fabbri

source("./code/multivariate_analysis/model.R")
source("./code/multivariate_analysis/preprocess.R")
source("./code/multivariate_analysis/dictionaries.R")

library(readr)
library(dplyr)
library(purrr)
library(ggraph)
library(tibble)
library(tidygraph)
library(parallel)
library(corpcor)
library(GeneNet)
library(Matrix)
library(qgraph)
library(gtsummary)
library(tidyselect)
library(data.table)
library(fastDummies)
library(sqldf)

#####################################################
##### Function performing main pipeline for GMM #####
#####################################################
main.pipeline.ggm <- function(params, save.data) {
  
  # Delete images in directory
  # dir.delete <- paste0("./results/ggm/", params$key.save, "/")
  # files.del <- dir(path = dir.delete, pattern = "*.png")
  # file.remove(file.path(dir.delete, files.del))
  # files.del <- dir(path = dir.delete, pattern = "*.pdf")
  # file.remove(file.path(dir.delete, files.del))
  # files.del <- dir(path = dir.delete, pattern = "*.csv")
  # file.remove(file.path(dir.delete, files.del))
  
  ## Parameters
  # Dictionary to create correct model.code
  codes <- params.to.model.code(model.type = "GMM")
  code.scaling.type <- codes$code.scaling.type
  code.adj.type <- codes$code.adj.type
  code.model.type <- codes$code.model.type
  #code.time.type
  code.omic.type <- codes$code.omic.type
  code.exposure.type <- codes$code.exposure.type
  
  # Parameters
  if (length(params$omic.type) == length(code.omic.type) - 1) {
    omic.types <- c("all")
  }
  stratification.groups <- names(codes$code.exposure.type) %>%
    .[1:length(.) - 1]
  stratification <- FALSE
  stratify.by.outcome <- FALSE
  stratification.outcome <- c("low", "normal", "high")
  time.point <- c(1, 2)
  package.corr <- params$package.corr
  method.na <- "kNN"
  perform.adj <- list(perform = TRUE, method = "residual")
  scale.new <- list(perform = TRUE, group = "none", 
                    method = "autoscale")
  plotting <- list(perform = FALSE, save = TRUE)
  validate <- FALSE
  is.hpc <- params$is.hpc
  
  # Create dataframe with combinations of parameters
  if (stratification == FALSE) {
    stratification.groups <- "all"
  }
  if (stratify.by.outcome == FALSE) {
    stratification.outcome <- "none"
  }
  all.params <- list(omic.type = omic.types, 
                     stratification.group = stratification.groups, 
                     stratification.outcome = stratification.outcome, 
                     adjustment.confounders = perform.adj$method, 
                     time.point = time.point)
  all.params <- expand.grid(all.params)
  print(all.params)
  
  ## Analyses
  # Iterate over combinations (rows) of parameters and perform analysis
  list.fitted.models <- list()
  for (row in 1:nrow(all.params)) {
    el <- list()
    for (param in names(all.params)) {
      el[[param]] <- all.params[row, param]
    }
    el$method.na   <- method.na
    el$scale.new   <- scale.new
    el$perform.adj <- perform.adj
    el$plotting    <- plotting
    el$stratification <- stratification
    el$stratify.by.outcome <- stratify.by.outcome
    el$boot <- params$boot
    
    # scaling . adjustment . correlation . time . omic . exposure
    el$model.code <- paste0("mod_", 
                            code.scaling.type[[el$scale.new$method %>% as.character()]], ".", 
                            code.adj.type[[el$perform.adj$method %>% as.character()]], ".", 
                            #code.model.type[[el$type.correlation %>% as.character()]], ".", 
                            el$time.point, ".", 
                            code.omic.type[[el$omic.type %>% as.character()]], ".", 
                            code.exposure.type[[el$stratification.group %>% as.character()]])
    #el$path.save <- path.save
    el$validate <- validate
    el$is.hpc <- is.hpc
    
    # Method/package to be used to compute the correlation matrix (i.e., the GGM)
    el$package.corr <- package.corr
    
    res <- perform.analysis(params = el, save.data = save.data)
    gc()
    
    list.fitted.models[[el$model.code]][["model"]] <- res$ggm
    list.fitted.models[[el$model.code]][["params"]] <- el
    list.fitted.models[[el$model.code]][["sample.size"]] <- res$sample.size
    list.fitted.models[[el$model.code]][["data"]] <- res$data
    
  } # End loop over parameters
  
  return(list.fitted.models)
} # End function main pipeline

#################################################################
##### Function to perform analysis given list of parameters #####
#################################################################
perform.analysis <- function(params, save.data) {
  # Load all data: -omics, exposures, metadata
  if (params$is.hpc == TRUE) {
    ## We are using data from the new HPC cluster
    
    # Filtered CG sites
    filtered.cgs <- readr::read_csv("data/cpgs_ewaffSVA_common.csv", 
                                    col_names = TRUE, 
                                    col_types = cols()) %>%
      as.list() %>% .$value
    
    path.data.methylome <- "../data/methylome/"
    name.methylome <- "methylome_panel_ComBatSlide_6cells_v4.csv"
    
    if (params$boot$perform) {
      # Load only selected CpG sites for bootstrapping
      selected.biomarkers <- read.csv("data/merged_biomarkers_ggm.csv") %>%
        tibble::as_tibble()
      cpg.sites <- selected.biomarkers %>%
        dplyr::filter(label == "methylome") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(name.ready = strsplit(name.ready, "_")[[1]][1])
      meth <- data.table::fread(paste0(path.data.methylome, name.methylome), 
                                select = c("SampleID", "HelixID", 
                                           cpg.sites$name.ready), 
                                nThread = 10) %>%
        dplyr::as_tibble()
    } else {
      meth <- data.table::fread(paste0(path.data.methylome, name.methylome), 
                                nThread = 10) %>%
        dplyr::as_tibble()
    }
    ############################################################################
    gc()
    meth <- meth %>%
      dplyr::filter(., grepl(ifelse(
        params$time.point == 1, "_1A", "_1B"
      ), SampleID)) %>%
      dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
      dplyr::mutate(SubjectID = substr(HelixID, 4, nchar(.))) %>%
      dplyr::select(c(SampleID, HelixID, tidyselect::starts_with("cg"))) %>%
      dplyr::mutate(dplyr::across(tidyselect:::where(is.factor), as.numeric)) %>%
      dplyr::rename_with(., ~ paste0(., "_me") %>% unname() %>% unlist(), 
        tidyselect::starts_with("cg"))

    # Filter CpG sites
    samples <- meth %>% dplyr::select(c(SampleID, HelixID))
    meth <- meth %>%
      dplyr::select(-c(SampleID, HelixID))
    if (!params$boot$perform) {
      meth <- meth %>%
      dplyr::select(tidyselect::all_of(filtered.cgs))
    }
    meth <- meth %>%
      minfi::logit2(.) %>%
      dplyr::bind_cols(samples)
    gc()
    
    path.meta <- "data/"
    path.data <- "../data/"
    
  } else {
    ## We are using data from the old cluster
    
    path.meta <- "data/"
    path.data <- "../data/"
  }

  exposome <- readRDS(file = paste0(path.data, "exposome_1", 
                                    ifelse(params$time.point == 1, 
                                            "A", "B"))) %>%
    .$data %>% tibble::as_tibble() %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID))

  prot <- readRDS(file = paste0(path.data, "proteome_1", 
                                ifelse(params$time.point == 1, 
                                        "A", "B")))
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                  ifelse(params$time.point == 1, 
                                          "A", "B")))
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                  ifelse(params$time.point == 1, 
                                          "A", "B")))

  list.omics.data <- list(prot$data, metabs$data, metabu$data)
  gc()
  if (params$is.hpc == TRUE) { list.omics.data <- append(list.omics.data, 
                                                         list(meth)) }
  omic <- list.omics.data %>%
    purrr::reduce(dplyr::inner_join, by = c("SampleID")) %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID))
  
  rm(list = c("prot", "metabs", "metabu"))
  if (params$is.hpc == TRUE) { rm("meth") }
  gc()
  
  season <- readr::read_csv(paste0(path.meta, "metadata_old.csv"), 
                            col_names = TRUE, col_types = cols()) %>%
    dplyr::select(c(SampleID, season, period)) %>%
    dplyr::filter(period == ifelse(params$time.point == 1, "A", "B")) %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::select(-c(period)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(HelixID = substr(SampleID, 1, nchar(SampleID) - 3)) %>%
    dplyr::select(-c(SampleID))
  
  metadata <- readr::read_csv(file = paste0(path.meta, "meta", 
                                            params$time.point, ".csv"), 
                              col_names = TRUE, col_types = cols()) %>%
    dplyr::select(-tidyselect::contains("HelixID")) %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID)) %>%
    dplyr::inner_join(season)
  
  # Filter out subjects not present in other time point and/or data types
  common.subjects <- readr::read_csv(paste0(path.meta, 
                                            paste0("common_samples_t", 
                                                   params$time.point, ".csv")), 
                                     col_types = cols()) %>%
    `colnames<-`(c("HelixID"))
  
  ######################### Bootstrapping #########################
  # Bootstrapping for subjects
  if (params$boot$perform) {
    ## Subjects from t1
    common.subjects1 <- readr::read_csv(paste0(path.meta, 
                                              paste0("common_samples_t", 1, ".csv")), 
                                        col_types = cols()) %>%
      `colnames<-`(c("HelixID"))
    ## Subjects from t2
    common.subjects2 <- readr::read_csv(paste0(path.meta, 
                                              paste0("common_samples_t", 2, ".csv")), 
                                        col_types = cols()) %>%
      `colnames<-`(c("HelixID"))
    ## Common subjects by time point
    common.subjects <- base::intersect(common.subjects1$HelixID, 
                                       common.subjects2$HelixID)
    
    ## Sample subjects with replacement
    seed <- params$boot$seed # Same by time point, different for each bootstrapping
    set.seed(seed = seed)
    common.subjects <- sample.int(n = length(common.subjects), replace = TRUE)
    common.subjects %>%
      write.csv(file = paste0(params$boot$path.save.subjects, 
                              "subjects_", seed, ".csv"), 
                quote = FALSE, row.names = FALSE)
  }
  ##############################################################################
  
  if (!params$boot$perform) {
    exposome <- exposome %>%
      dplyr::filter(HelixID %in% common.subjects$HelixID) %>%
      dplyr::select(-c("SampleID"))
    omic <- omic %>%
      dplyr::filter(HelixID %in% common.subjects$HelixID) %>%
      dplyr::select(-c("SampleID"))
    metadata <- metadata %>%
      dplyr::filter(HelixID %in% common.subjects$HelixID) %>%
      dplyr::rename(zBMI = hs_zbmi_theano) %>%
      dplyr::select(-c("SampleID"))
  } else {
    exposome <- exposome[common.subjects, ] %>%
      dplyr::select(-c("SampleID"))
    omic <- omic[common.subjects, ] %>%
      dplyr::select(-c("SampleID"))
    metadata <- metadata[common.subjects, ] %>%
      dplyr::rename(zBMI = hs_zbmi_theano) %>%
      dplyr::select(-c("SampleID"))
  }
  gc()
  
  dat <- list(
    exposures = exposome, 
    omics     = omic, 
    metadata  = metadata
  )
  if (!params$boot$perform & save.data) {
    base::saveRDS(object = dat$exposures, 
                  file = "../data/dist_vars/exposures_raw")
    base::saveRDS(object = dat$omics, 
                  file = "../data/dist_vars/omics_raw")
    base::saveRDS(object = dat$metadata, 
                  file = "../data/dist_vars/metadata_raw")
  }
  
  # Run main analysis
  ggm.fitted <- fit.ggm(data = dat, params = params, 
                        save.data = save.data)
  rm(dat)
  gc()
  
  return(list(
    ggm = ggm.fitted$corr.res, 
    sample.size = ggm.fitted$sample.size, 
    data = ggm.fitted$data
  ))
} # End function perform analysis

#######################################
##### Function to fit GGM to data #####
#######################################
fit.ggm <- function(data, params, save.data) {
  
  ## Prepare data
  # Covariates for adjusting models
  covariates <- list.covariates()
  
  # Stratification by exposure family/group
  if (params$stratification) {
    stop(call. = TRUE)
    data$exposures <- data$exposures %>%
      dplyr::select(., dplyr::contains(c("HelixID", 
                                         params$stratification.group %>%
                                           as.character())))
  }
  # Stratification by outcome
  if (params$stratify.by.outcome) {
    stop(call. = TRUE)
  }
  
  # Adjust for covariates
  if (params$perform.adj$perform == TRUE) {
    if (params$perform.adj$method == "variables") {
      all.covariates <- list.covariates()
      metadata.tmp <- fastDummies::dummy_cols(.data = data$metadata[, 
                                                                    all.covariates$covariates.char], 
                                              select_columns = all.covariates$covariates.char, 
                                              remove_first_dummy = FALSE, 
                                              remove_selected_columns = TRUE)
      data$metadata <- cbind(metadata.tmp, data$metadata[, 
                                                         all.covariates$covariates.num])
      colnames(data$metadata) <- lapply(colnames(data$metadata), function(x) {
        paste0(x, "_co")
      })
      
      data$exposures <- cbind(data$exposures, data$metadata)
    } else {
      omics <- perform.adj(data = data$omics %>% dplyr::select(-c(HelixID)), 
                           metadata = data$metadata, 
                           covariates = covariates, 
                           method = params$perform.adj$method)
      data$omics <- base::cbind(data$omics[, c("HelixID")], 
                                omics)
      colnames(data$omics)[1] <- "HelixID"
    }
  }
  gc()
  
  # Transform datasets
  data <- transform.data(data = data, scale.new = params$scale.new, 
                         is.adjusted = params$perform.adj$perform)
  
  ## Fit model
  data <- base::cbind(data$exposures, data$omics)
  
  ######################### Bootstrapping #########################
  # Filter selected biomarkers from merged network
  nrow.old <- nrow(data)
  if (params$boot$perform) {
    selected.biomarkers <- read.csv("data/merged_biomarkers_ggm.csv") %>%
      .$name.ready
    data <- as.data.frame(data) %>%
      dplyr::select(dplyr::any_of(selected.biomarkers))
    
    # Case when not using Methylation data
    if (params$is.hpc == FALSE) {
      selected.biomarkers <- selected.biomarkers[!grepl("_me", 
                                                        selected.biomarkers)]
    }
    
    if (dim(data)[2] != length(selected.biomarkers)) { stop(call. = TRUE) }
    if (nrow(data) != nrow.old) { stop(call. = TRUE) }
  }
  gc()
  ##############################################################################
  
  if (!params$boot$perform & save.data) {
    base::saveRDS(object = data, file = "../data/dist_vars/data_matrix_processed")
  }
  dim.data <- dim(data)[1]
  cat("\n##################################################\n")
  cat(paste0("Dimension of dataset: ", 
             dim(data)[1], "x", dim(data)[2], ".\n"))
  num.chemicals <- data %>%
    dplyr::select(tidyselect::contains("_e"))
  cat(paste0("Number of chemicals: ", dim(num.chemicals)[2], "."))
  cat("\n##################################################\n")
  
  if (params$package.corr == "corr.corr") {
    stop()
    corr.res <- correlation::correlation(data = data, 
                                         method = params$type.correlation %>% as.character(), 
                                         p_adjust = params$type.padjust %>% as.character(), 
                                         partial = params$partial.corr)
    
  } else if (params$package.corr == "corpcor") {
    # (Shrinkage) Partial correlation coefficients
    corr.res <- corpcor::pcor.shrink(x = data)
    gc()
    
  } else if (params$package.corr == "huge.glasso") {
    corr.res <- bootnet::estimateNetwork(data = data, default = "huge", 
                                         tuning = 0.1, graphType = "pcor", 
                                         verbose = TRUE,  datatype = "normal", 
                                         nCores = 6)
    
  } else if (params$package.corr == "vanilla.huge") {
    corr.res <- huge::huge(x = data %>% as.matrix(), 
                           method = "glasso", nlambda = 50)
  }
  
  gc()
  ret <- list(
    corr.res = corr.res, 
    sample.size = dim.data, 
    data = data
  )
  return(ret)
} # End function to fit GGMs

#################################################################
##### Function to process results of GGM fitting procedures #####
#################################################################
process.ggms <- function(list.ggms, active, 
                         filter.mixed.interactions, is.directed, 
                         params, boot = FALSE, save.data) {
  
  processed.res <- list()
  for (cor.res in names(list.ggms)) {
    
    # Different packages return different results w/ different
    # APIs. Here we retrieve the same information from all of them
    if (active == "corpcor") {
      # Number of correlation coefficients computed
      len <- attributes(list.ggms[[cor.res]]$model)$dim[1]
      labels.rows <- attributes(list.ggms[[cor.res]]$model)$dimnames[[1]]
      labels.columns <- attributes(list.ggms[[cor.res]]$model)$dimnames[[2]]
      # Estimated shrinkage coefficient
      lambda <- attributes(list.ggms[[cor.res]]$model)$lambda
      
      # "Raw" correlation matrix
      pcorr.mat <- list.ggms[[cor.res]]$model
      
    } else if (active == "huge.glasso") {
      # Dataframe with correlation coefficients
      df <- list.ggms[[cor.res]]$model$graph %>%
        `diag<-`(1.0)
      colnames(df) <- list.ggms[[cor.res]]$model$labels
      df <- df %>% tibble::as_tibble() %>%
        dplyr::mutate(row.names = colnames(.)) %>%
        tibble::column_to_rownames(., var = "row.names")
      
      # "Raw" correlation matrix
      pcorr.mat <- list.ggms[[cor.res]]$model$graph %>%
        `diag<-`(1.0)
      
      # Labels to replace ID nodes
      labels.rows <- list.ggms[[cor.res]]$model$labels
      labels.columns <- list.ggms[[cor.res]]$model$labels
      
    } else if (active == "relevance") {
       pcorr.mat <- list.ggms[[cor.res]]$model
       
       # In this case, pcorr.mat is *not* a symmetric matrix, 
       # which is required by the function network.test.edges
       pcorr.mat.sym <- diag(dim(pcorr.mat)[1] + dim(pcorr.mat)[2])
       
       labels.rows <- rownames(pcorr.mat)
       labels.columns <- colnames(pcorr.mat)
       rownames(pcorr.mat.sym) <- c(labels.rows, labels.columns)
       colnames(pcorr.mat.sym) <- c(labels.rows, labels.columns)
       labels.rows <- rownames(pcorr.mat.sym)
       labels.columns <- colnames(pcorr.mat.sym)
       
       pcorr.mat.sym[row.names(pcorr.mat), colnames(pcorr.mat)] <- pcorr.mat
       pcorr.mat.sym <- Matrix::forceSymmetric(pcorr.mat.sym, uplo = "U") %>%
         as.matrix()
       if (!isSymmetric(pcorr.mat.sym)) { stop(call. = TRUE) }
       
       # Helper function to check whether the symmetric matrix is 
       # same as original data. Very slow if large matrix!!!
       check.symmetric.mat <- function(sym, orig) {
         for (r in row.names(orig)) {
           for (c in colnames(orig)) {
             val.sym <- sym[r, c]
             val.orig <- orig[r, c]
             if (val.sym != val.orig) { return(FALSE) }
           }
         }
         
         return(TRUE)
       }
       
       #checked <- check.symmetric.mat(pcorr.mat.sym, pcorr.mat)
       #if (checked == FALSE) { stop(call. = TRUE) }
       
       pcorr.mat <- pcorr.mat.sym
      
    } else if (active == "vanilla.huge") {
      # Select optimal path
      selected <- huge::huge.select(est = list.ggms[[cor.res]]$model, 
                                    criterion = "stars")
      
      # "Raw" correlation matrix
      pcorr.mat <- selected$icov[[selected$opt.index]]
      
      # Labels to replace ID nodes
      labels.columns <- colnames(list.ggms[[cor.res]]$model$data)
      labels.rows <- labels.columns
      
      # Dataframe with correlation coefficients
      df <- selected$icov[[selected$opt.index]]
      colnames(df) <- labels.columns
      df <- df %>% tibble::as_tibble() %>%
        dplyr::mutate(row.names = colnames(.)) %>%
        tibble::column_to_rownames(., var = "row.names")
    }
    
    # Store correlation matrix
    if (boot == FALSE & save.data) {
      base::saveRDS(object = pcorr.mat, 
                    file = paste0("../data/correlations/pcorr_mat_raw_", 
                                  cor.res))
    }
    
    ## Compute p-values for partial correlation coefficients
    df.significance <- GeneNet::network.test.edges(r.mat = pcorr.mat, 
                                                   fdr = TRUE, cutoff.method = "fndr", 
                                                   direct = FALSE, plot = FALSE)
    old.dim <- dim(pcorr.mat)[1]
    rm(pcorr.mat)
    gc()
    if (dim(df.significance)[1] != old.dim * (old.dim - 1) / 2) { stop(call. = TRUE) }
    
    # Replace ID nodes w/ actual labels
    if ((length(labels.rows) != length(labels.columns)) | (!identical(labels.rows, 
                                                                      labels.columns))) {
      stop(call. = TRUE)
      }
    df.significance <- df.significance %>%
      dplyr::rowwise() %>%
      dplyr::mutate(node1 = labels.rows[node1]) %>%
      dplyr::mutate(node2 = labels.columns[node2])
    gc()
    if (boot == FALSE & save.data) {
      base::saveRDS(object = df.significance, 
                    file = paste0("../data/correlations/network_raw_", 
                                  cor.res))
    }
    
    # Eventually filter network based on either probability or q-value
    if (params$method.ggm == "none") {
      df.filtered <- df.significance %>%
        # Reorder columns
        dplyr::select(., c(node1, node2, pcor, pval, qval, prob))
    } else {
      df.filtered <- df.significance %>%
        GeneNet::extract.network(network.all = ., 
                                 method.ggm = params$method.ggm, 
                                 cutoff.ggm = params$cutoff.ggm) %>%
        # Reorder columns
        dplyr::select(., c(node1, node2, pcor, pval, qval, prob))
    }
    rm(df.significance)
    gc()
    
    # Tidy network for visualization
    df.filtered <- df.filtered %>%
      dplyr::filter(pcor != 0) %>%
      tidygraph::as_tbl_graph(., directed = is.directed) %>%
      tidygraph::activate(., what = "nodes") %>%
      # Community detection
      dplyr::mutate(community = as.factor(tidygraph::group_infomap())) %>%
      # Centrality
      dplyr::mutate(degree = tidygraph::centrality_degree()) %>%
      tidygraph::activate(., what = "edges") %>%
      # Direction edges based on sign of correlation values
      dplyr::mutate(direction = as.factor(sign(pcor)))
    
    # Add type for visualizing difference in nodes
    df.filtered <- df.filtered %>%
      tidygraph::activate(., what = "nodes") %>%
      dplyr::mutate(type = name %>% sapply(., function(x) {
        strsplit(x, "_")[[1]][2]
      }))
    type.nodes <- df.filtered %>% dplyr::select(type) %>% tibble::as_tibble() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(type2 = ifelse(
        type %in% c("ms", "mu", "p", "me", "e"), type, "co"))
    df.filtered <- df.filtered %>%
      dplyr::mutate(type = type.nodes$type2)
    
    # Store results for current model
    #processed.res[[cor.res]][["corr"]] <- df # Correlation matrix
    processed.res[[cor.res]][["net"]]  <- df.filtered # Actual network w/ p-values
    gc()
    
    if (boot == FALSE & save.data) {
      base::saveRDS(object = df.filtered, 
                    file = paste0("../data/correlations/network_filtered_", 
                                  cor.res))
    }
    
    # Save final network
    gg <- ggraph::ggraph(processed.res[[cor.res]][["net"]], 
                         layout = "kk") +
      ggraph::geom_edge_link(colour = "darkgrey", alpha = 0.75) +
      ggraph::geom_node_point(mapping = aes(shape = type, size = degree, 
                                            color = type)) +
      ggraph::geom_node_text(mapping = aes(label = name), repel = TRUE)
    if (boot == FALSE) {
      ggplot2::ggsave(paste0(params$path.save.plots, cor.res, ".pdf"), 
                      dpi = 720/2, 
                      width = 20, height = 12)
    }
  } # End loop fitted models
  
  return(processed.res)
} # End function post-process correlation matrices (i.e., GGMs)

##################################################
##### Function to compute network properties #####
##################################################
net.properties <- function(ggms, type.networks, 
                           path.save, code.save = "", 
                           is.merged, how.to.join = NULL) {
  
  # Iterate over fitted models
  properties.all <- list()
  df.all <- list()
  for (n in names(ggms)) {
    l <- ggms[[n]]$net
    
    # code GGMs: scaling . adjustment . time . omic . exposure group
    if (type.networks == "ggm") {
      code <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
        as.list() %>% .[[5]]
    } else if (type.networks == "relevance") {
      code <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
        as.list() %>% .[[1]] %>% strsplit(., "") %>% .[[1]] %>% .[6]
    }
    
    node.attributes <- l %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
    edge.attributes <- l %>% tidygraph::activate(., what = "edges") %>%
      tidygraph::as_tibble()
    
    dim.old <- dim(edge.attributes)[1]
    df <- dplyr::inner_join(node.attributes, edge.attributes, 
                            by = c("idx" = "from"))
    df <- dplyr::inner_join(node.attributes, df, 
                            by = c("idx" = "to"))
    
    if (dim(df)[1] != dim(edge.attributes)[1]) { stop(call. = TRUE) }
    
    if (is.merged == FALSE) {
      colnames(df) <- c("from", "community.x", "degree.x", "type.x", "idx.x", 
                        "to", "community.y", "degree.y", "type.y", "idx.y", 
                        "pcor", "pval", "qval", "prob", "direction")
      df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                pcor, pval, qval, direction))
    } else {
      if (how.to.join == "inner") {
        colnames(df) <- c("from", "type.x", "group.x", "idx.x", 
                          "to", "type.y", "group.y", "idx.y", 
                          "type.x.2", "type.y.2", "pcor.x", "pval.x", 
                          "qval.x", "direction", "pcor.y", "pval.y", "qval.y")
        df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                  pcor.x, pcor.y, 
                                  pval.x, pval.y, 
                                  qval.x, qval.y, direction))
      } else {
        colnames(df) <- c("from", "type.x", "idx.x", "to", "type.y", 
                          "idx.y", "type.x.2", "type.y.2", 
                          "pcor", "pval", "qval", "direction")
        df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                  pcor, pval, qval, direction))
      }
    }
    
    df$code <- as.character(n)
    if (dim.old != dim(df)[1]) {stop(call. = TRUE)}
    
    df.all <- append(df.all, list(df))
    
  } # End loop over networks
  
  # Create summary of properties of interest
  df <- df.all %>%
    purrr::reduce(rbind)
  properties <- df %>%
    dplyr::select(-c(from, to)) %>%
    gtsummary::tbl_summary(by = "code") %>%
    gtsummary::as_gt() %>%
    gt::gtsave(filename = paste0(path.save, "props_nets", code.save, ".pdf"), 
               vwidth = 1500)
  
  return(properties)
} # End function compute properties networks

############################################################################
##### Function to merge networks based on parameter (e.g., time point) #####
############################################################################
merge.networks <- function(ggms, exposure.group, omic.type = "none", 
                           how.to.join, type.networks, path.save, 
                           boot = FALSE, save.data) {
  
  # Iterate over fitted models
  df.all <- list()
  for (n in names(ggms)) {
    l <- ggms[[n]]$net
    
    # code: scaling . adjustment . time . omic . exposure group
    if (type.networks == "ggm") {
      code <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
        as.list() %>% .[[5]]
    } else if (type.networks == "relevance") {
      code <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
        as.list() %>% .[[1]] %>% strsplit(., "") %>% .[[1]] %>% .[6]
    }
    if (code != exposure.group) next
    
    node.attributes <- l %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
    edge.attributes <- l %>% tidygraph::activate(., what = "edges") %>%
      tidygraph::as_tibble()
    
    dim.old <- dim(edge.attributes)[1]
    df <- dplyr::inner_join(node.attributes, edge.attributes, 
                            by = c("idx" = "from"))
    df <- dplyr::inner_join(node.attributes, df, 
                            by = c("idx" = "to"))
    if (dim(df)[1] != dim(edge.attributes)[1]) { stop(call. = TRUE) }
    colnames(df) <- c("from", "community.x", "degree.x", "type.x", "idx.x", 
                      "to", "community.y", "degree.y", "type.y", "idx.y", 
                      "pcor", "pval", "qval", "prob", "direction")
    df <- dplyr::select(df, c(from, to, type.x, type.y, 
                              pcor, pval, qval, direction))
    if (dim.old != dim(df)[1]) {stop(call. = TRUE)}
    
    df.all <- append(df.all, list(df))
  } # End loop over GGMs
  
  # Check percentage of inter-connections for later comparison 
  # with merged network
  df1Info <- summarize.net(df.all[[1]] %>% tidygraph::as_tbl_graph())
  df2Info <- summarize.net(df.all[[2]] %>% tidygraph::as_tbl_graph())
  
  df.info <- sqldf::sqldf(x = "
                            SELECT *
                            FROM df1Info
                            JOIN df2Info
                              ON df1Info.`type.a` = df2Info.`type.a` AND df1Info.`type.b` = df2Info.`type.b` OR
                                df1Info.`type.a` = df2Info.`type.b` AND df1Info.`type.b` = df2Info.`type.a`
                          ")
  colnames(df.info) <- c("type.a", "type.b", 
                         "n.1", "p.1", 
                         "type.a2", "type.b2", 
                         "n.2", "p.2")
  df.info <- df.info %>%
    dplyr::select(-c(type.a2, type.b2))
  
  # Merge dfs for the 2 time points to keep only shared edges
  # In this case, the `.x` or `.y` in the merged df refer to 
  # the 2 time points and not the two nodes
  if (how.to.join == "inner") {
    df.1 <- dplyr::inner_join(df.all[[1]], df.all[[2]], 
                            by = c("from" = "from", "to" = "to", 
                                   "direction"))
    df.2 <- dplyr::inner_join(df.all[[1]], df.all[[2]], 
                            by = c("from" = "to", "to" = "from", 
                                   "direction"))
    df <- dplyr::bind_rows(df.1, df.2)
    
    df.check <- df %>% dplyr::filter(type.x.x != type.x.y)
    if (!(identical(df.check$type.x.x, df.check$type.y.y) & identical(df.check$type.y.x, df.check$type.x.y))) {
      stop(call. = TRUE)
    }
    
    df$type.x.y <- NULL
    df$type.y.y <- NULL
    df <- df %>%
      dplyr::rename(type.x = type.x.x, type.y = type.y.x)
    
  } else if (how.to.join == "anti.x") {
    stop(call. = TRUE)
    # Connections in t1 *not* present in t2
    df <- dplyr::anti_join(df.all[[1]], df.all[[2]], 
                           by = c("from", "to", "direction", 
                                  "type.x", "type.y"))
  } else if (how.to.join == "anti.y") {
    stop(call. = TRUE)
    # Connections in t2 *not* present in t1
    df <- dplyr::anti_join(df.all[[2]], df.all[[1]], 
                           by = c("from", "to", "direction", 
                                  "type.x", "type.y"))
  } else {
    stop(call. = TRUE)
  }
  gc()
  
  # Store merged network
  if (boot == FALSE & save.data) {
    base::saveRDS(object = df, file = "../data/correlations/merged_net")
  }
  
  # Plot merged networks
  net <- df %>% tidygraph::as_tbl_graph(directed = TRUE)
  net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.))) %>%
    dplyr::rowwise() %>%
    # Add labels
    dplyr::mutate(label = strsplit(name, "_")[[1]][2])
  
  net <- net %>%
    tidygraph::to_undirected() %>%
    tidygraph::activate(., what = "nodes") %>%
    tidygraph::mutate(label = factor(net.nodes$label)) %>%
    tidygraph::mutate(group = tidygraph::group_components())
  
  # Summarize connections in final network
  if (omic.type != "none") { exposure.group <- paste0("o", omic.type, 
                                                      "e", exposure.group)}
  
  info.merged <- summarize.net(net)
  info <- sqldf::sqldf(x = "
                        SELECT *
                        FROM `info.merged` AS m
                        JOIN `df.info` as b
                          ON m.`type.a` = b.`type.a` AND m.`type.b` = b.`type.b` OR
                          m.`type.a` = b.`type.b` AND m.`type.b` = b.`type.a`
                       ")
  colnames(info) <- c("type.a", "type.b", 
                      "n", "p", 
                      "type.ao", "type.bo", 
                      "n.1", "p.1", "n.2", "p.2")
  info <- info %>%
    dplyr::select(-c(type.ao, type.bo)) %>%
    dplyr::arrange(type.a, type.b, desc(n), desc(p)) %>% print()
  
  pie.chart <- info %>%
    dplyr::select(-c(p, p.1, p.2)) %>%
    tidyr::gather("time", "count", -type.a, -type.b) %>%
    dplyr::mutate(time = replace(time, time == "n", "merged")) %>%
    dplyr::mutate(time = replace(time, time == "n.1", "t1")) %>%
    dplyr::mutate(time = replace(time, time == "n.2", "t2")) %>%
    dplyr::group_by(time) %>%
    dplyr::mutate(count.cumsum = count / sum(count)) %>%
    dplyr::mutate(count.cumsum.cumsum = cumsum(count.cumsum)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., mapping = aes(x = count.cumsum, 
                                     y = as.factor(time), 
                                     fill = interaction(type.a, type.b, 
                                                        sep = "-"))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar() +
    ggplot2::theme_bw() +
    ggplot2::geom_text(aes(label = count), 
                       position = position_stack(vjust = 0.6)) +
    ggplot2::theme(panel.border = element_blank(), panel.grid = element_blank(), 
                   axis.ticks = element_blank(), 
                   axis.text.x = element_blank(), axis.text.y = element_text(size = 25), 
                   axis.title = element_blank(), 
                   legend.text = element_text(size = 25), 
                   legend.title = element_text(size = 25)) +
    ggplot2::labs(fill = "edge type") +
    ggplot2::scale_fill_brewer(palette = "Spectral")
  # ggplot2::ggsave(paste0("output/posters/pptox22_tex/images/piechart_edges.pdf"), 
  #                 dpi = 720/2, 
  #                 width = 20, height = 12)
  
  # info %>%
  #   gt::gt() %>%
  #   gt::gtsave(., file = "output/posters/pptox22_tex/dat/summary_edges.tex")
  
  # Recode labels to long format and remove coding from variables' names
  net <- net %>%
    dplyr::mutate(label = dplyr::recode(
      label, 
      e = "exposure", 
      ms = "serum metabolome", 
      mu = "urinary metabolome", 
      p = "proteome", 
      me = "methylome"
    )) %>%
    dplyr::mutate(name = ifelse(
      grepl("log.", name), 
      stringr::str_replace(name, "log.", ""), name
    )) %>%
    tidygraph::activate(., what = "nodes") %>%
    tidygraph::mutate(name = ifelse(
      label == "exposure", toupper(name), name
    )) %>%
    dplyr::mutate(name = stringr::str_replace(name, "_.*", ""))
  
  gg <- ggraph::ggraph(net, 
                       layout = "kk") +
    ggraph::geom_edge_link(colour = "darkgrey", alpha = 0.75) +
    ggraph::geom_node_point(mapping = aes(shape = label, 
                                          color = label), 
                            size = 4) +
    ggraph::geom_node_text(mapping = aes(label = name), 
                           repel = TRUE, size = 5) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   axis.title = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), 
                   axis.ticks = ggplot2::element_blank(), 
                   legend.title = element_text(size = 25), 
                   legend.text = element_text(size = 22))
  if (boot == FALSE) {
    ggplot2::ggsave(paste0(path.save, "merged_", exposure.group, ".pdf"), 
                    dpi = 720/2, 
                    width = 20, height = 12)
  }
  
  # Plot distribution of correlation coefficients from 2 time points 
  # and the relative change for each edge/association
  if (how.to.join == "inner") {
    vec.pcor <- c("pcor.x", "pcor.y")
  } else { vec.pcor <- c("pcor") }
  
  dist.coeffs <- net %>%
    tidygraph::activate(., what = "edges") %>%
    tidygraph::as_tibble() %>% dplyr::select(dplyr::all_of(vec.pcor)) %>%
    tidyr::gather() %>%
    ggplot2::ggplot(aes(x = value, fill = key)) +
    ggplot2::geom_histogram(alpha = 0.6, position = "identity") +
    ggplot2::labs(title = exposure.group)
  if (boot == FALSE) {
    ggplot2::ggsave(paste0(path.save, "distCoeffs_", exposure.group, ".pdf"), 
                    dpi = 720/2, 
                    width = 20, height = 12)
  }
  
  edges <- net %>%
    tidygraph::activate(., what = "edges") %>%
    tidygraph::as_tibble()
  dim.old <- dim(edges)[1]
  nodes <- net %>%
    tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>%
    dplyr::mutate(id = seq.int(nrow(.)))
  
  edges <- dplyr::inner_join(edges, nodes, by = c("from" = "id"))
  edges <- dplyr::inner_join(edges, nodes, by = c("to" = "id"))
  if (dim(edges)[1] != dim.old) { stop(call. = TRUE) }
  
  ret <- NULL
  if (omic.type == "none") { exposure.group <- n }
  ret[[exposure.group]] <- list("net" = net)
  return(ret)
} # End function merge networks

################################################################################
##### Helper function to create file with molecule information for network #####
################################################################################
join.mol.info <- function(nodes, info, nodes.col, info.col, 
                          path.save, file.name, idx) {
  
  ret <- dplyr::inner_join(nodes, info, 
                           by = setNames(info.col, nodes.col))
  if (dim(ret)[1] > 0) {
    readr::write_csv(ret, file = paste0(path.save, "cc_", idx, "_", 
                                        file.name, "_mols.csv"), 
                     col_names = TRUE)
  }
  
}

###############################################################
##### Function to find and plot exposure-omic connections #####
###############################################################
find.edges.exposures <- function(path.df, cc.idx, path.save) {
  # Input: path to .csv file with results of connected component
  
  all.exposures <- dict.exposure.groups() %>%
    unname() %>% toupper()
  all.edges <- readr::read_csv(path.df, col_types = cols())
  colnames(all.edges) <- c("node_a", "node_b")
  
  df.filtered <- all.edges %>%
    dplyr::filter((node_a %in% all.exposures & !(node_b %in% all.exposures)) |
                    (node_b %in% all.exposures & !(node_a %in% all.exposures)))
  
  # Create graph for direct connections only
  label.mols <- readr::read_csv("data/labels_mols.csv", 
                                col_names = TRUE, col_types = cols()) %>%
    dplyr::mutate(node.clean = ifelse(
      label == "e", toupper(node.clean), node.clean
    ))
  label.mols$node <- NULL
  df.filtered <- dplyr::inner_join(df.filtered, label.mols, 
                                   by = c("node_a" = "node.clean"))
  df.filtered <- dplyr::inner_join(df.filtered, label.mols, 
                                  by = c("node_b" = "node.clean"))
  colnames(df.filtered) <- c("from", "to", "type.x", "type.y")
  net <- tidygraph::as_tbl_graph(df.filtered, directed = FALSE)
  
  # Add labels nodes
  net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
  net.edges <- net %>% tidygraph::activate(., what = "edges") %>%
    tidygraph::as_tibble() %>%
    dplyr::select(c(from, to, type.x, type.y))
  
  for (idx in 1:nrow(net.nodes)) {
    label <- net.edges[net.edges$from == idx, ]$type.x
    if (length(label) == 0) {
      label <- net.edges[net.edges$to == idx, ]$type.y
    }
    
    if (length(label) > 1 & length(unique(label)) == 1) {
      net.nodes[idx, "label"] <- label[1]
    } else if (length(label) == 1) {
      net.nodes[idx, "label"] <- label
    } else if (length(label) > 1 & length(unique(label)) > 1) {
      stop(call. = TRUE)
    }
  }
  
  net <- net %>%
    tidygraph::activate(., what = "nodes") %>%
    tidygraph::mutate(label = factor(net.nodes$label)) %>%
    dplyr::mutate(label = dplyr::recode(
      label, 
      e = "exposure", 
      ms = "serum metabolome", 
      mu = "urinary metabolome", 
      p = "proteome", 
      me = "methylome"
    ))
  
  gg <- ggraph::ggraph(net, 
                       layout = "kk") +
    ggraph::geom_edge_link(colour = "darkgrey", alpha = 0.75) +
    ggraph::geom_node_point(mapping = aes(shape = label, 
                                          color = label), 
                            size = 4) +
    ggraph::geom_node_text(mapping = aes(label = name), 
                           repel = TRUE, size = 5) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   axis.title = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), 
                   axis.ticks = ggplot2::element_blank(), 
                   legend.title = element_text(size = 25), 
                   legend.text = element_text(size = 22))
  ggplot2::ggsave(paste0(path.save, "focus_cc_", cc.idx, ".pdf"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
  
  return(net)
}

####################################################
##### Function to analyse connected components #####
####################################################
analyse.cc <- function(net, path.save) {
  
  # Load raw data to retrieve information on molecules
  path.data <- "../data/"
  
  prot <- readRDS(file = paste0(path.data, "proteome_1A")) %>% .$feature.data
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1A")) %>% .$feature.data %>%
    dplyr::mutate(CHEBI = as.character(CHEBI)) %>%
    dplyr::mutate(CHEBI = ifelse(grepl("/", CHEBI), 
                                 strsplit(CHEBI, "/")[[1]][1], 
                                 CHEBI)) %>%
    dplyr::select(c(Rvar, var, Class, CHEBI))
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1A")) %>% .$feature.data %>%
    dplyr::mutate(CHEBI = as.character(CHEBI)) %>%
    dplyr::mutate(CHEBI = ifelse(grepl("/", CHEBI), 
                                 strsplit(CHEBI, "/")[[1]][1], 
                                 CHEBI)) %>%
    dplyr::select(c(Rvar, var, CHEBI, KEGG))
  
  grouped.net <- net %>%
    tidygraph::group_by(group)
  num.cc <- grouped.net %>% dplyr::select(group) %>% tidygraph::as_tibble() %>%
    max()
  # Iterate over found connected components
  for (idx in 1:num.cc) {
    # Filter one connected component
    tmp <- grouped.net %>%
      tidygraph::filter(group == idx)
    
    ## Create df with molecules of current connected component for further analysis
    mols <- tmp %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble()
    # If connected component w/ exposures, add label to file name
    #if ("exposure" %in% mols$label) { idx <- paste0(idx, "_e") }
    
    mols <- mols %>%
      dplyr::filter(!label %in% c("e")) %>% dplyr::select(name) %>%
      dplyr::mutate(mol = name %>% sapply(., function(x) {
        strsplit(x, split = "_")[[1]][1]
      })) %>% dplyr::select(mol)
    
    if (dim(mols)[1] > 0) {
      join.mol.info(nodes = mols, info = prot, 
                    nodes.col = "mol", info.col = "Prot_ID", 
                    path.save = path.save, file.name = "prot", idx = idx)
      join.mol.info(nodes = mols, info = metabs, 
                    nodes.col = "mol", info.col = "Rvar", 
                    path.save = path.save, file.name = "metabs", idx = idx)
      join.mol.info(nodes = mols, info = metabu, 
                    nodes.col = "mol", info.col = "Rvar", 
                    path.save = path.save, file.name = "metabu", idx = idx)
    }
    
    # "Visualize" table of merged network
    net.nodes <- tmp %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
    net.edges <- tmp %>% tidygraph::activate(., what = "edges") %>%
      tidygraph::as_tibble()
    
    info <- dplyr::inner_join(net.edges, net.nodes, 
                              by = c("from" = "idx"))
    info <- dplyr::inner_join(info, net.nodes, 
                              by = c("to" = "idx"))
    info <- info %>%
      dplyr::select(c(name.x, name.y, pcor.x, pcor.y, 
                      qval.x, qval.y, pval.x, pval.y)) %>%
      dplyr::mutate(name.x = name.x %>% sapply(., function(x) {
        strsplit(x, split = "_")[[1]][1] %>%
          ifelse(grepl("log.", .), 
                 strsplit(., split = "log.")[[1]][2], 
                 .)
      })) %>%
      dplyr::mutate(name.y = name.y %>% sapply(., function(x) {
        strsplit(x, split = "_")[[1]][1] %>%
          ifelse(grepl("log.", .), 
                 strsplit(., split = "log.")[[1]][2], 
                 .)
      }))
      readr::write_csv(x = info %>%
                         dplyr::select(c(name.x, name.y)), 
                       file = paste0(path.save, "cc_", idx, ".csv"), 
                       col_names = TRUE)
    
    # Plot
    pcor.means <- tmp %>% tidygraph::activate(., what = "edges") %>%
      tidygraph::as_tibble() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pcor.median = median(pcor.x, pcor.y))
    
    phenols <- c("MEPA", "ETPA", "PRPA", "BUPA", "BPA", "TCS", "OXBE")
    phthalates <- c("MEP", "MIBP", "MNBP", "MBZP", "MEHP", "MEHHP", 
                    "MEOHP", "MECPP", "OHMINP", "OXOMINP")
    op.pesticides <- c("DMP", "DMTP", "DMDTP", "DEP", "DETP")
    
    tmp %>%
      tidygraph::activate(., what = "nodes") %>%
      tidygraph::mutate(name = ifelse(
        label == "exposure", toupper(name), name
      )) %>%
      tidygraph::mutate(class = ifelse(
        name %in% phenols, "phenols", 
        ifelse(
          name %in% phthalates, "phthalates", 
          ifelse(
            name %in% op.pesticides, "OP pesticides", ""
          )
        )
      )) %>%
      tidygraph::activate(what = "edges") %>%
      tidygraph::mutate(pcor.median = pcor.means$pcor.median) %>%
      ggraph::ggraph(graph = ., 
                     layout = "kk") +
      ggforce::geom_mark_hull(mapping = aes(x = x, y = y, 
                                            fill = class, 
                                            filter = class != ""), 
                              concavity = 1) +
      ggplot2::scale_fill_brewer(palette = "Dark2") +
      ggraph::geom_edge_link(colour = "darkgrey", 
                             mapping = aes(width = pcor.median), 
                             show.legend = FALSE) +
      ggraph::geom_node_point(mapping = aes(shape = label, 
                                            color = label), 
                              size = 8) +
      ggraph::geom_node_text(mapping = aes(label = name), 
                             repel = TRUE, nudge_x = 0.06, nudge_y = 0.06, 
                             size = 10, 
                             show.legend = FALSE) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.line = ggplot2::element_blank(), 
                     axis.title = ggplot2::element_blank(), 
                     axis.text = ggplot2::element_blank(), 
                     axis.ticks = ggplot2::element_blank(), 
                     legend.text = element_text(size = 25), 
                     legend.title = element_text(size = 25)) +
      ggplot2::labs(fill = "")
    ggplot2::ggsave(paste0(path.save, "cc_", idx, ".pdf"), 
                    dpi = 720/2, 
                    width = 20, height = 12)
  }
  
}

########################################################################################################
##### Function to extract molecule names from connected components for enrichment/pathway analysis #####
########################################################################################################
prepare.molecules <- function() {
  molecules.chebi <- readr::read_tsv("data/compounds.tsv", col_types = cols())
  path.dat <- "results/ggm/"
  
  # Iterate over files with molecules in each connected component
  list.ccs <- list.files("results/ggm/", pattern = "*mols\\.csv") %>%
    tibble::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cc = as.numeric(strsplit(value, "_") %>% .[[1]] %>% .[2])) %>%
    dplyr::mutate(type = strsplit(value, "_") %>% .[[1]] %>% .[3]) %>%
    dplyr::filter(type != "prot") %>%
    dplyr::arrange(cc)
  
  # Iterate over connected components
  for (cc.idx in 1:max(list.ccs$cc)) {
    df <- list.ccs %>%
      dplyr::filter(cc == cc.idx)
    if (dim(df)[1] == 0) { next }
    
    # Merge molecules in connected component with data from CHEBI DB
    dat <- readr::read_csv(paste0(path.dat, df$value[1]), 
                           col_types = cols()) %>%
      dplyr::mutate(CHEBI = paste0("CHEBI:", CHEBI))
    if (dim(df)[1] > 1) {
      dat.2 <- readr::read_csv(paste0(path.dat, df$value[2]), 
                               col_types = cols()) %>%
        dplyr::mutate(CHEBI = paste0("CHEBI:", CHEBI))
      
      dat <- dplyr::bind_rows(dat, dat.2) %>%
        dplyr::inner_join(molecules.chebi, by = c("CHEBI" = "CHEBI_ACCESSION"))
    } else {
      dat <- dat %>%
        dplyr::inner_join(molecules.chebi, by = c("CHEBI" = "CHEBI_ACCESSION"))
    }
    
    cat(dat$NAME, sep = "\n")
    cat("\n##########\n")
  }
  
}

###########################################
##### Functions to summarize networks #####
###########################################
count.edge.by.node <- function(g) {
  df <- data.frame(igraph::get.edgelist(g))
  res <- aggregate(rep(1, nrow(df)), by = as.list(df), FUN = sum)
  names(res) <- c("type.a", "type.b", "n")
  res <- res %>%
    dplyr::mutate(p = round(n / sum(n), 2) * 100) %>%
    dplyr::arrange(type.a, type.b, desc(n), desc(p))
  
  return(res)
}

summarize.net <- function(net) {
  
  info <- net %>%
    tidygraph::activate(., what = "edges") %>%
    tidygraph::as_tibble() %>%
    dplyr::select(c(type.x, type.y))
  
  igraph.obj <- igraph::graph.data.frame(info, directed = FALSE)
  info <- count.edge.by.node(igraph.obj)
  
  return(info)
}
