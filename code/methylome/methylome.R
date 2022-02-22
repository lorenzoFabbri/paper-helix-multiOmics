# Functions to load, pre-process and clean Methylome data
# Author: Lorenzo Fabbri

library(tidyverse)
library(methylumi)
library(parallel)
library(sfsmisc)
library(ewaff)

source("./code/multivariate_analysis/dictionaries.R")

###############################################################
##### Function to load Methylome and filter by time point #####
###############################################################
load.methylome <- function(filter.time) {
  path.data.methylome <- "../data/methylome/"
  #name.methylome <- "methylome_panel_ComBatSlide_6cells_notfitr_v4.csv"
  name.methylome <- "methylome_panel_ComBatSlide_6cells_v4.csv"
  
  # Load data
  meth <- readr::read_csv(file = paste0(path.data.methylome, name.methylome), 
                          num_threads = 8, show_col_types = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::mutate(SubjectID = substr(HelixID, 4, nchar(.))) %>%
    dplyr::select(c(SampleID, HelixID, tidyselect::starts_with("cg"))) %>%
    dplyr::mutate(dplyr::across(tidyselect:::where(is.factor), as.numeric)) %>%
    dplyr::rename_with(., ~ paste0(., "_me") %>% unname() %>% unlist(), 
                       tidyselect::starts_with("cg"))
  
  # Eventually filter by time point
  if (filter.time$perform == TRUE) {
    meth <- meth %>%
      dplyr::filter(., grepl(ifelse(
        filter.time$time.point == 1, "_1A", "_1B"
      ), SampleID))
  }
  
  gc()
  return(meth)
}

####################################################################
##### Filter Methylome to exclude near-zero variance CpG sites #####
####################################################################
filter.by.variability <- function(dat, cutoff) {
  betas <- dat %>%
    dplyr::select(tidyselect::contains("_me")) %>%
    as.matrix() %>% t()
  gc()
  betas <- new("MethyLumiSet", betas = betas)
  
  # Filter out low-variability sites
  betas.filtered <- methylumi::varFilter(betas, var.func = IQR, 
                                         var.cutoff = cutoff)$eset %>%
    Biobase::exprs() %>% as.matrix() %>% t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(SampleID = dat$SampleID) %>%
    dplyr::mutate(HelixID = dat$HelixID)
  
  if (dim(betas.filtered)[1] != dim(dat)[1]) { stop(call. = TRUE) }
  
  cat(paste0("Keeping top ", (1 - cutoff) * 100, "% features.\n"))
  cat(paste0("Original data: ", dim(dat)[2], ".\n"))
  cat(paste0("Filtered data: ", dim(betas.filtered)[2], "."), "\n")
  
  rm(list = c("betas"))
  gc()
  return(betas.filtered)
}

##################################################
##### Filter Methylome based on results EWAS #####
##################################################
filter.by.ewas <- function(dat, time.point) {
  ## Formula: beta_i ~ exposure_j + covariates, for all i's and j's
  
  path.data <- "../data/"
  path.meta <- "data/"
  
  # Load exposures
  exposome <- readRDS(file = paste0(path.data, "exposome_1", ifelse(
    time.point == 1, "A", "B"
  ))) %>%
    .$data %>% tibble::as_tibble() %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID))
  
  # Load covariates
  season <- readr::read_csv(paste0(path.meta, "metadata_old.csv"), 
                            col_names = TRUE) %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::select(c(SampleID, season, period)) %>%
    dplyr::filter(period == ifelse(
      time.point == 1, "A", "B"
    )) %>%
    dplyr::select(-c(period))
  
  metadata <- readr::read_csv(file = paste0(path.meta, "meta", 
                                            time.point, ".csv"), 
                              col_names = TRUE) %>%
    dplyr::rename(zBMI = hs_zbmi_theano) %>%
    dplyr::select(-tidyselect::contains("HelixID")) %>%
    dplyr::inner_join(season, by = "SampleID")
  
  # Filter only common subjects between Exposome and Methylome
  common.subjects <- intersect(intersect(dat$SampleID, exposome$SampleID), 
                               metadata$SampleID)
  
  cat("\nFiltering data by SampleID...\n")
  dat <- dat %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  exposome <- exposome %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  metadata <- metadata %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  cat("\n... Done.\n")
  if ((dim(dat)[1] != dim(exposome)[1]) | (dim(dat)[1] != dim(metadata)[1])) {
    #if (dim(dat)[1] != dim(exposome)[1]) {
    stop(call. = TRUE)
  }
  
  # Prepare datasets for analyses
  dat <- dat %>% dplyr::select(tidyselect::contains("_me"))
  exposome <- exposome %>% dplyr::select(-c(SampleID, HelixID))
  exposome <- scale(exposome, center = TRUE, scale = TRUE)
  metadata <- metadata %>% dplyr::select(dplyr::all_of(
    list.covariates() %>% unlist() %>% unname()
  )) %>%
    dplyr::mutate(dplyr::across(where(is.character), as_factor))
  
  X <- cbind(exposome, metadata)
  #X <- exposome
  gc()
  
  # Perform EWAS (iterate over betas and for each beta over each exposure)
  ## Helper function to fit one exposure and one CpG site
  fit.single.ewas <- function(beta.val, x) {
    beta.val <- tibble::as_tibble(beta.val)
    x <- tibble::as_tibble(x)
    
    # Robust lm using MASS
    fit <- MASS::rlm(formula = as.formula(paste(colnames(beta.val), "~ .")), 
                     data = cbind(beta.val, x))
    summ.fit <- data.frame(summary(fit)$coefficients)
    
    # Computing p-value for the exposure using robust Wald's test
    summ.fit$p.val <- survey::regTermTest(fit, test.terms = colnames(x), 
                                          #df = summary(fit)$df[2], 
                                          null = NULL, df = Inf, 
                                          method = "Wald")$p[1]
    summ.fit$p.val2 <- survey::regTermTest(fit, test.terms = colnames(x), 
                                           df = summary(fit)$df[2], 
                                           null = NULL, 
                                           method = "Wald")$p[1]
    summ.fit$n.obs <- stats::nobs(fit)
    
    return(summ.fit[2, ] %>% tibble::rownames_to_column(var = "exposure"))
  }
  
  ewas.res <- list()
  for (exposure in colnames(exposome)) {
    cat(paste0("Fitting models for: ", exposure, ".\n"))
    x <- X %>% dplyr::select(dplyr::all_of(c(exposure, colnames(metadata))))
    #x <- X %>% dplyr::select(dplyr::all_of(exposure))
    
    # Iterate over the betas
    summ.fit <- parallel::mclapply(X = dat, # beta values
                                   FUN = fit.single.ewas, 
                                   x = x, 
                                   mc.cores = 19)
    #mc.cores = parallel::detectCores() - 1)
    ewas.res <- append(ewas.res, summ.fit)
    
    gc()
  }
  
  return(dplyr::bind_rows(ewas.res, .id = "cpg"))
}

#################################################
##### Filter Methylome based on reliability #####
#################################################
filter.by.reliability <- function(dat, cutoff.rel) {
  
  reliability.results <- openxlsx::read.xlsx("../data/methylome/Sugden_MethylationReliability_Data_S1.xlsx", 
                                             startRow = 3) %>%
    dplyr::filter(Reliability > cutoff.rel)
  
  # Methylation sites form HELIX data
  meth <- dat %>% dplyr::select(tidyselect::contains("_me")) %>%
    colnames(.) %>%
    lapply(., function(x) {
      strsplit(x, "_") %>% unlist() %>% .[1]
    }) %>% unlist()
  
  # Find intersection
  reliable.sites.helix <- base::intersect(reliability.results$Illumina.Probe.ID, 
                                          meth)
  
  return(reliable.sites.helix)
}

##############################################
##### Main function to run EWAS w/ EWAFF #####
##############################################
methylome.ewaff <- function(meth, time.point, perform.adj) {
  options(mc.cores = 25)
  path.data <- "../data/"
  path.meta <- "data/"
  
  ewas.res <- list()
  ewas.res[[time.point]] <- list()
  
  # Filter time point
  meth.time <- meth %>%
    dplyr::filter(., grepl(ifelse(
      time.point == 1, "_1A", "_1B"
    ), SampleID))
  
  # Load exposures
  exposome <- readRDS(file = paste0(path.data, "exposome_1", ifelse(
    time.point == 1, "A", "B"
  ))) %>%
    .$data %>% tibble::as_tibble() %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID))
  
  # Load covariates
  season <- readr::read_csv(paste0(path.meta, "metadata_old.csv"), 
                            col_names = TRUE, show_col_types = FALSE) %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::select(c(SampleID, season, period)) %>%
    dplyr::filter(period == ifelse(
      time.point == 1, "A", "B"
    )) %>%
    dplyr::select(-c(period))
  metadata <- readr::read_csv(file = paste0(path.meta, "meta", 
                                            time.point, ".csv"), 
                              col_names = TRUE, show_col_types = FALSE) %>%
    dplyr::rename(zBMI = hs_zbmi_theano) %>%
    dplyr::select(-tidyselect::contains("HelixID")) %>%
    dplyr::inner_join(season, by = "SampleID")
  
  # Filter only common subjects between Exposome and Methylome
  common.subjects <- intersect(intersect(meth.time$SampleID, exposome$SampleID), 
                               metadata$SampleID)
  meth.time <- meth.time %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  exposome <- exposome %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  metadata <- metadata %>% dplyr::filter(SampleID %in% common.subjects) %>%
    dplyr::arrange(SampleID)
  if ((dim(meth.time)[1] != dim(exposome)[1]) | (dim(meth.time)[1] != dim(metadata)[1])) {
    stop(call. = TRUE)
  }
  
  # Filter by reliability score (0.6=good, 0.75=excellent)
  reliable.sites <- filter.by.reliability(dat = meth.time, 
                                          cutoff.rel = 0.6) %>%
    paste0(., "_me")
  meth.time <- meth.time %>%
    dplyr::select(dplyr::all_of(c("SampleID", reliable.sites)))
  
  # Prepare datasets for analyses
  meth.time <- meth.time %>% dplyr::select(tidyselect::contains("_me")) %>%
    t()
  gc()
  exposome <- exposome %>% dplyr::select(-c(SampleID, HelixID))
  exposome <- robustHD::robStandardize(exposome)
  metadata <- metadata %>% dplyr::select(dplyr::all_of(
    list.covariates() %>% unlist() %>% unname()
  )) %>%
    dplyr::mutate(dplyr::across(where(is.character), as_factor))
  
  # Fit models
  if (perform.adj == TRUE) {
    X <- cbind(exposome, metadata)
  } else {
    X <- exposome %>% tibble::as_tibble()
  }
  for (exposure in colnames(exposome)) {
    cat(paste0("Fitting models for: ", exposure, ".\n"))
    
    # Winsorize exposure
    #X[[exposure]] <- robustHD::winsorize(X[[exposure]], standardized = TRUE)
    
    if (perform.adj == TRUE) {
      x <- X %>% dplyr::select(dplyr::all_of(c(exposure, colnames(metadata))))
    } else {
      x <- X %>% dplyr::select(exposure)
    }
    x <- data.frame(x)
    gc()
    
    # Handle putative outliers
    methylation.no.outliers <- ewaff::ewaff.handle.outliers(methylation = meth.time, 
                                                            method = "winsorize")[[1]]
    colnames(methylation.no.outliers) <- paste0("s", 
                                                1:dim(methylation.no.outliers)[2])
    cpg.sites <- rownames(methylation.no.outliers)
    methylation.no.outliers <- as.data.frame(methylation.no.outliers)
    rownames(methylation.no.outliers) <- cpg.sites
    
    ############################################################################
    # Otherwise complains when computing SVs (e.g. SVA)
    methylation.no.outliers <- as.matrix(methylation.no.outliers)
    gc()
    ############################################################################
    
    if (perform.adj == TRUE) {
      formula.tmp <- paste("methylation~", exposure, "+", 
                           paste(colnames(metadata), collapse = "+"))
    } else {
      formula.tmp <- paste("methylation~", exposure)
    }
    
    res.tmp <- ewaff::ewaff.sites(formula = formula.tmp, 
                                  variable.of.interest = exposure, 
                                  methylation = methylation.no.outliers, 
                                  data = x, 
                                  method = "rlm", 
                                  generate.confounders = "smartsva")
    
    ret <- res.tmp$table %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "cpg") %>%
      dplyr::as_tibble()
    ret$exposure <- exposure
    
    # Store tibble w/ results fitted model
    ret %>%
      readr::write_csv(., paste0("results/ewaff/t", time.point, 
                                 "_", exposure, 
                                 ".csv"), 
                       col_names = TRUE)
    
    ewas.res[[time.point]] <- append(ewas.res[[time.point]], list(ret))
    gc()
  } # End loop exposures
  
  rm(list = c("meth.time", "exposome", "season", "metadata"))
  gc()
  
  rm(meth)
  gc()
  
  return(ewas.res)
} # End function EWAS w/ EWAFF

###########################################################################
##### Function to load results EWAFF (by time point and by chemicals) #####
###########################################################################
load.res.ewaff <- function(path.res) {
  
  # Helper function to load all files given list of file names
  .load.files <- function(path, list.file.names) {
    ret <- list()
    for (el in list.file.names) {
      tmp <- vroom::vroom(paste0(path, el), col_types = cols(), 
                          num_threads = 6)
      ret <- append(ret, list(tmp))
    }
    
    ret <- dplyr::bind_rows(ret)
    return(ret)
  }
  
  t1 <- list.files(path = path.res, pattern = "t1_")
  t2 <- list.files(path = path.res, pattern = "t2_")
  
  # Data frames containing all pair-wise associations, by time point
  t1 <- .load.files(path.res, t1)
  t2 <- .load.files(path.res, t2)
  
  return(list(t1, t2))
}

#############################################
##### Function to process results EWAFF #####
#############################################
process.res.ewaff <- function(df1, df2, threshold.fdr, key.save) {
  
  # Adjust p-values by chemical using FDR instead of Bonferroni
  df1 <- df1 %>%
    dplyr::group_by(exposure) %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup()
  df2 <- df2 %>%
    dplyr::group_by(exposure) %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup()
  
  # Find intersection of CpGs w/ filtering on FDR-adjusted p-values
  cpgs.1 <- df1 %>%
    dplyr::filter(fdr <= threshold.fdr) %>%
    dplyr::select(cpg)
  cpgs.2 <- df2 %>%
    dplyr::filter(fdr <= threshold.fdr) %>%
    dplyr::select(cpg)
  cpgs.filtered <- base::intersect(
    cpgs.1 %>% c() %>% .[["cpg"]], 
    cpgs.2 %>% c() %>% .[["cpg"]]
  ) %>% unique()
  
  # Filter dataframes based on common CpGs
  df1.filtered <- df1 %>%
    dplyr::filter(cpg %in% cpgs.filtered)
  df2.filtered <- df2 %>%
    dplyr::filter(cpg %in% cpgs.filtered)
  
  cpgs.filtered %>%
    tibble::as_tibble() %>%
    readr::write_csv(., paste0("data/cpgs_ewaff", toupper(key.save), "_common.csv"), 
                     col_names = TRUE)
  
  return(list(
    "cpgs.filtered" = cpgs.filtered, 
    "df1.filtered" = df1.filtered, 
    "df2.filtered" = df2.filtered
  ))
}

#############################################
##### Main function to run all analyses #####
#############################################
methylome.main <- function(meth, key.save, perform.ewas) {
  
  # Filter sites by "reliability"
  # Reference: https://www.sciencedirect.com/science/article/pii/S2666389920300143
  # With cutoff=0.6, we select >="good" sites
  reliable.sites <- filter.by.reliability(dat = meth, cutoff.rel = 0.6) %>%
    paste0(., "_me")
  meth.reliable <- meth %>%
    dplyr::select(dplyr::all_of(c("SampleID", reliable.sites)))
  cat(paste0("Number of samples: ", dim(meth.reliable)[1], ".\n"))
  cat(paste0("Dimension of Methylome after reliability filtering: ", dim(meth.reliable)[2], ".\n"))
  cat(paste0("Number of removed probes after reliability filtering: ", dim(meth)[2] - dim(meth.reliable)[2], ".\n"))
  
  meth.reliable %>%
    readr::write_csv(., paste0("./data/cpgs_rel_", 
                               key.save, ".csv"), 
                     col_names = TRUE)
  gc()
  
  if (perform.ewas) {
    for (time.point in c(1, 2)) {
      # Filter time point
      meth.time <- meth.reliable %>%
        dplyr::filter(., grepl(ifelse(
          time.point == 1, "_1A", "_1B"
        ), SampleID))
      
      # Filter sites by EWAS w/ exposures
      meth.ewas <- filter.by.ewas(dat = meth.time, time.point = time.point)
      cat(paste0("Number of removed probes after EWAS: ", dim(meth.time)[2] - dim(meth.ewas)[2], ".\n"))
      
      meth.ewas %>%
        readr::write_csv(., paste0("./data/cpgs_t", time.point, "_ewas_", 
                                   key.save, ".csv"), 
                         col_names = TRUE)
      
      rm(list = c("meth.time", "meth.ewas"))
      gc()
    } # End loop time points
  } # EWAS end
  
}

#################################################
##### Function to post-process results EWAS #####
#################################################
postprocess.ewas <- function(path.df.ewas, threshold.padj, time.point) {
  
  # Code for Meff
  # exposome <- readRDS(file = paste0(path.data, "exposome_1", ifelse(
  #   time.point == 1, "A", "B"
  # ))) %>%
  #   .$data %>% tibble::as_tibble() %>%
  #   dplyr::mutate(HelixID = from.sample.to.helix(SampleID)) %>%
  #   dplyr::select(-c(SampleID, HelixID))
  # exposome <- scale(exposome, center = TRUE, scale = TRUE)
  # cormat <- exposome %>%
  #   stats::cor(use = "pairwise.complete.obs")
  # lambdas <- eigen(cormat)$values
  # Meff <- ncol(cormat) - sum((lambdas > 1) * (lambdas - 1))
  ##############################################################################
  
  if (typeof(path.df.ewas) == "character") {
    df <- readr::read_csv(path.df.ewas)
  } else {
    df <- path.df.ewas
  }
  
  df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cpg.id = sapply(cpg, function(x) {
      strsplit(x, "_")[[1]][1]
    })) %>%
    dplyr::group_by(exposure) %>%
    dplyr::mutate(p.adj = p.adjust(p.val, method = "fdr")) %>%
    dplyr::ungroup() #%>%
  #dplyr::mutate(p.adj.meff = pmin(1, Meff * p.val))
  #dplyr::filter(p.adj <= threshold.padj)
  colnames(df) <- c("cpg", "exposure", "beta", "sde", "t", "p", "cpg.id", 
                    "fdr", "fdr.meff")
  
  unique.cpgs <- unique(df$cpg.id)
  return(df)
}

##### Helper function to merge CpG sites from two time points
merge.cpgs.ewas <- function(df1, df2) {
  ret <- unique(base::intersect(df1$cpg, df2$cpg))
  
  return(ret)
}
