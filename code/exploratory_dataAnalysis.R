# Author: Lorenzo Fabbri
# Functions to perform Exploratory Data Analysis

library(dplyr)
library(GGally)
library(ggplot2)
library(ggcorrplot)
library(heatmaply)
library(ComplexHeatmap)
library(tidyHeatmap)
suppressMessages(library(ggh4x, quietly = TRUE))
library(reshape2)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(knitr)

##### HELPER FUNCTIONS #####

source("./code/prepare_datasets.R")

# Count number of missing values by column
na_byCol <- function(dat) {
  as.data.frame(colSums(is.na(get(dat))))
}

load.data <- function(time.period) {
  dat <- prepare.datasets(time.period = time.period)
  return(dat)
}

# Correlation within exposures
correlation.variables <- function(variables) {
  round.num <- 2
  corr.data <- round(cor(variables), round.num)
  return(corr.data)
}

# Correlation of exposures between time points
correlation.variables.time <- function(dat1, dat2) {
  
  # Filter subjects present in both time point measurements and 
  # remove missing values
  dat1 <- tidyr::drop_na(dat1)
  dat2 <- tidyr::drop_na(dat2)
  common <- dplyr::inner_join(dat1, dat2, "SubjectID")
  
  tmp.data1 <- dat1 %>%
    dplyr::filter(SubjectID %in% common$SubjectID) %>%
    dplyr::select(where(is.numeric))
  tmp.data2 <- dat2 %>%
    dplyr::filter(SubjectID %in% common$SubjectID) %>%
    dplyr::select(where(is.numeric))
  
  round.num <- 2
  corr.data <- round(cor(tmp.data1, tmp.data2), round.num)
  
  return(corr.data)
}
correlation.variables.time.run <- function(type.data) {
  if (type.data == 'exposome') {
    names.data <- c('exp_metab_blood', 'exp_metabun', 'exp_proteome')
    col <- c("OP Pesticides" = "Red", 
             "Phenols" = "Blue", 
             "Phthalates" = "Green")
  } else if (type.data == 'omics') {
    names.data <- c('omic_metab_blood', 'omic_metabun', 'omic_proteome')
    col <- c("omic_metab_blood" = "Red", 
             "omic_metabun" = "Blue", 
             "omic_proteome" = "Green")
  }
  
  # Correlation analysis (between time points)
  for (datum in names.data) {
    # This is necessary only to obtain the right variable labels
    tmp.data <- get(datum) %>%
      dplyr::select(where(is.numeric)) %>%
      tidyr::drop_na()
    
    corr.data <- correlation.variables.time(get(datum), 
                                            get(paste0(datum, '2')))
    ha <- ComplexHeatmap::rowAnnotation(foo = anno_text(var_label(tmp.data), 
                                                        gp = gpar(col = col)))
    plt <- ComplexHeatmap::Heatmap(corr.data, 
                                   cluster_rows = T, #left_annotation = ha, 
                                   show_row_dend = F, show_column_dend = F, show_row_names = F, 
                                   width = 500, height = 200, 
                                   use_raster = T, 
                                   # row_names_max_width = max_text_width(text = unlist(var_label(tmp.data)), 
                                   #                                      gp = gpar(fontsize = 2)), 
                                   name = paste('Heatmap (between time points): ', datum))
    print(plt)
  }
}

eda.exposome <- function(time.period) {
  exposomes <- c('exp_metab_blood', 'exp_metabun', 'exp_proteome')
  if (time.period == 2) {
    exposomes <- paste(exposomes, '2', sep = '')
  }
  col <- c("OP Pesticides" = "Red", 
           "Phenols" = "Blue", 
           "Phthalates" = "Green")
  
  # Check for NaN
  cat('Checking for NaN...\n')
  print(lapply(exposomes, na_byCol))
  
  # Check distribution exposures
  cat('Checking distributions...\n')
  tmp.data <- get(exposomes[2]) %>%
    dplyr::select(where(is.numeric)) %>%
    gather() %>%
    ggplot(aes(value)) +
      geom_density() +
      suppressMessages(stat_theodensity(color = "blue"), classes = "message") +
      facet_wrap(~key, scales = "free") +
      ggtitle(paste('Distribution variables in Exposome. Time point: ', time.period))
  print(tmp.data)
  
  # Correlation analysis exposures (within time point)
  cat('Checking for correlations...\n')
  for (exposome in exposomes[1]) {
    # Drop unnecessary columns
    tmp.data <- get(exposome) %>%
      dplyr::select(where(is.numeric))
    
    corr.data <- correlation.variables(tmp.data)
    ha <- ComplexHeatmap::rowAnnotation(foo = anno_text(var_label(tmp.data), 
                                                        gp = gpar(col = col)))
    c.heatmap <- ComplexHeatmap::Heatmap(corr.data, 
                                         cluster_rows = T, left_annotation = ha, 
                                         show_column_dend = F, show_row_names = F, 
                                         width = 100, height = 100, 
                                         use_raster = T, 
                                         row_names_max_width = max_text_width(text = unlist(var_label(tmp.data)), 
                                                                              gp = gpar(fontsize = 2)), 
                                         name = paste('Heatmap: ', exposome, 
                                                      '. Time point: ', time.period))
    print(c.heatmap)
  }
  
  # PCA
  cat('Performing PCA...\n')
  res.pca <- FactoMineR::PCA(tmp.data, 
                             scale. = T, graph = F)
  factoextra::fviz_pca_biplot(res.pca, 
                              label = "var", 
                              col.var = unlist(var_label(tmp.data)), 
                              title = paste('PCA of Exposome. Time point: ', time.period))
}

eda.omics <- function(time.period) {
  names.data <- c('omic_metab_blood', 'omic_metabun', 'omic_proteome')
  if (time.period == 2) {
    names.data <- paste(names.data, '2', sep = '')
  }
  
  # Check for NaN
  cat('Checking for NaN...\n')
  print(lapply(names.data, na_byCol))
  
  cat('Dropping rows containing missing values...\n')
  for (datum in names.data) {
    tmp <- get(datum) %>%
      tidyr::drop_na()
    assign(datum, tmp)
  }
  
  # Check distributions
  cat('Checking distributions...\n')
  for (idx in 1:3) {
    tmp.data <- get(names.data[idx]) %>%
      dplyr::select(where(is.numeric)) %>%
      gather() %>%
      ggplot(aes(value)) +
      geom_density() +
      suppressMessages(stat_theodensity(color = "blue"), classes = "message") +
      facet_wrap(~key, scales = "free") +
      ggtitle(paste('Distribution variables in: ', names.data[idx], '. Time point: ', time.period))
    
    # The figures are too big to be displayed inside RStudio
    path.plot <- paste0("/home/lorenzo/Downloads/", names.data[idx], "_dist", ".png")
    ggsave(path.plot, width = 20, height = 20)
    #print(tmp.data)
  }
  
  # Correlation analysis (within time point)
  cat('Checking for correlations...\n')
  for (datum in names.data) {
    # Drop unnecessary columns
    tmp.data <- get(datum) %>%
      dplyr::select(where(is.numeric))
    
    corr.data <- correlation.variables(tmp.data)
    path.plot <- paste0("/home/lorenzo/Downloads/", datum, "_corr", ".png")
    png(path.plot, width = 16, height = 16, units = "in", 
        res = 600)
    c.heatmap <- ComplexHeatmap::Heatmap(corr.data, 
                                         cluster_rows = T, show_column_dend = F, show_row_names = F, 
                                         #width = 100, height = 100, 
                                         use_raster = T, 
                                         name = paste('Heatmap: ', datum, 
                                                      '. Time point: ', time.period))
    print(c.heatmap)
    dev.off()
  }
  
  # PCA
  cat('Performing PCA...\n')
  for (datum in names.data) {
    tmp.data <- get(datum) %>%
      dplyr::select(where(is.numeric))
    
    res.pca <- FactoMineR::PCA(tmp.data, 
                               scale. = T, graph = F)
    plt <- factoextra::fviz_pca_biplot(res.pca, 
                                       label = "var", 
                                       title = paste('PCA of Omics: ', datum, 
                                                     '. Time point: ', time.period))
    print(plt)
  }
}
