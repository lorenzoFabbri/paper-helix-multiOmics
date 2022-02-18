# Author: Lorenzo Fabbri
# Function to merge/split datasets (-Omics, exposures, metadata)

library(readxl)
library(tidyverse)
library(janitor)
library(labelled)
library(magrittr)

prepare.datasets <- function(PATH = '/run/user/1000/gvfs/smb-share:server=fs01.isglobal.lan,share=data/ws_helix/HELIX_analyses/HELIX_Panel_Nonpersistent_multi-omics/', 
                             time.period) {
  ##### HELPER FUNCTIONS #####
  
  # Function to split dataframes with results ExWAS
  # into -Omics and Exposome
  split_data <- function(filename, sheet, time.period) {
    results <- readxl::read_xlsx(paste0(PATH, filename), 
                                 sheet = sheet)
    exposomeExternal <- unique(results$exposure)
    omics <- unique(results$omic)
    
    exposomeData <- rg %>%
      dplyr::select(c(HelixID, SubjectID, exposomeExternal))
    omicsData <- rg %>%
      dplyr::select(c(HelixID, SubjectID, omics))
    
    exposomeData <- change.type(exposomeData, "SubjectID", as.character)
    omicsData <- change.type(omicsData, "SubjectID", as.character)

    # Create dictionary of unique exposures
    dict <- create_dict(filename, sheet+1)
    # Merge information
    for (idx_col in colnames(exposomeData)[-c(1, 2)]) {
      label <- dict[dict$exposure == idx_col, 'subgroup.exposure']
      labelled::var_label(exposomeData[[idx_col]]) <- label[[1]]
    }
    
    out <- list(exposomeData, omicsData, exposomeExternal, omics)
    return(out)
  }
  
  # Function to create dictionary of variables (e.g. exposure family)
  create_dict <- function(filename, sheet) {
    cat("\t Creating dictionary... \n")
    results <- readxl::read_xlsx(paste0(PATH, filename), 
                                 sheet = sheet, .name_repair = "universal")
    
    # To simplify things, I'm gonna rename some columns...
    names(results) <- tolower(names(results))
    idxs_sameName <- which(names(results) == 'exposure')
    colnames(results)[idxs_sameName[1]] <- 'exposure.short'
    
    # Select only family, short and long names
    if (!('subgroup.exposure' %in% colnames(results))) {
      # I should not do this but I don't know whether I can 
      # modify the xlsx files
      idx <- which(names(results) == 'exp.group')
      colnames(results)[idx] <- 'subgroup.exposure'
    }
    results <- results %>%
      dplyr::select(c(subgroup.exposure, exposure.short, exposure)) %>%
      dplyr::distinct(exposure, .keep_all = T)
    
    return(results)
  }
  
  # Function to change type of given column
  change.type <- function(data, col, new.type) {
    data[[col]] %<>% new.type
    return(data)
  }
  
  ##### LOAD DATA #####
  
  # Load dataset containing all variables
  rg <- get(load(paste0(PATH, 'data/RG_sem', time.period, 'Sept2019.rdata')))
  
  # For some reason the same variables have different names
  # in the two datasets
  if (time.period == 2) {
    rg <- rg %>%
      dplyr::rename(., SubjectID = SubjectID.x)
  }
  
  ## Load datasets with ExWAS results to create dictionary variables
  for (type in c('metabun', 'metab_blood', 'proteome')) {
    cat(paste('Processing', type, sep=" "))
    cat('\n')
    
    tmp <- split_data(paste0('results/', type, '.xlsx'), 1, time.period)
    name1 <- paste('exp_', type, sep="")
    name2 <- paste('omic_', type, sep="")
    
    assign(name1, tmp[[1]])
    assign(name2, tmp[[2]])
    
    rg <- rg %>%
      dplyr::select(-dplyr::one_of(tmp[[4]]))
  }
  
  rg <- rg %>%
    dplyr::select(-dplyr::one_of(tmp[[3]]))
  
  # Check columns Exposomes are the same after ordering
  janitor::compare_df_cols_same(exp_metab_blood, exp_metabun, exp_proteome)
  
  # After filtering out the exposures and -Omics (columns),
  # there still are several non-Metadata columns left
  possible_exposures_left <- rg %>%
    dplyr::select(dplyr::starts_with(c('h_', 'hs_', 'hcp_', 'log.'))) %>%
    colnames()
  possible_metadata <- rg %>%
    dplyr::select(-dplyr::all_of(possible_exposures_left))
  
  ##### DATA PRE-PROCESSING #####
  
  ## Remove rows (i.e., subjects) with missing values
  omic_metabun <- omic_metabun %>%
    tidyr::drop_na()
  
  ## Change type columns containing "cdesc"
  exp_metab_blood$hcp_dmp_cdesc   <- as.factor(exp_metab_blood$hcp_dmp_cdesc)
  exp_metab_blood$hcp_dmdtp_cdesc <- as.factor(exp_metab_blood$hcp_dmdtp_cdesc)
  exp_metabun$hcp_dmp_cdesc       <- as.factor(exp_metabun$hcp_dmp_cdesc)
  exp_metabun$hcp_dmdtp_cdesc     <- as.factor(exp_metabun$hcp_dmdtp_cdesc)
  exp_proteome$hcp_dmp_cdesc      <- as.factor(exp_proteome$hcp_dmp_cdesc)
  exp_proteome$hcp_dmdtp_cdesc    <- as.factor(exp_proteome$hcp_dmdtp_cdesc)
  
  return(list(exp_metab_blood, exp_metabun, exp_proteome, 
              omic_metab_blood, omic_metabun, omic_proteome, 
              possible_metadata, possible_exposures_left))
}
