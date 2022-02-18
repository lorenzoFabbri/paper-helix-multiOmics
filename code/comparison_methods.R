# Author: Lorenzo Fabbri
# Functions to compare different regression methods

library(readr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggrepel)
library(rsample)
library(broom)
library(purrr)
library(glmnet)
library(fastDummies)
library(pls)
library(ggforce)

source("code/multivariate_analysis/preprocess.R")

# Helper function to perform 2-stage residual-outcome regression
perform.2sr <- function(outcome, covariates) {
  res.model <- lm(data.matrix(outcome) ~ data.matrix(covariates))
  predicted <- predict(res.model, covariates, type = "response")
  
  # Deflation step
  deflated <- outcome - predicted
  
  return(deflated)
}

# Function to prepare and pre-process data
generate.clean.data <- function(params) {
  time.point <- ifelse(params$time.point == 1, "", 2)
  
  exposome <- get(paste0("exp_", params$omic.type, time.point))
  omic     <- get(paste0("omic_", params$omic.type, time.point))
  metadata <- get(paste0("metadata", time.point))
  metadata_all <- readr::read_csv("./data/metadata_all.csv", col_names = TRUE)
  metadata <- metadata_all %>%
    dplyr::filter(period == ifelse(time.point == 1, "A", "B")) %>%
    dplyr::inner_join(metadata, by = "HelixID") %>%
    dplyr::mutate_if(sapply(., is.character), as.factor)
  
  dat <- list(exposures = exposome, omics = omic, metadata = metadata)
  
  # Pre-processing
  dat.preprocessed <- preprocess.data(data = dat, 
                                      scale.new = params$scale.new, 
                                      is.adjusted = params$adjust$perform, 
                                      validate = params$validate)
  
  # Adjusting for confounders, if any
  if (params$adjust$perform == TRUE) {
    if (params$adjust$method == "add.covariates") {
      covariates.names <- unlist(params$adjust$list, 
                                 recursive = TRUE, use.names = FALSE)
      formula.covariates <- paste(covariates.names, collapse = "+")
    } else if (params$adjust$method == "residual") {
      covariates.names <- unlist(params$adjust$list, 
                                 recursive = TRUE, use.names = FALSE) %>%
        list() %>% .[[1]]
      dat.preprocessed$omics <- perform.2sr(outcome = dat.preprocessed$omics, 
                                            covariates = dat.preprocessed$metadata[, covariates.names])
      formula.covariates <- ""
    }
  }
  
  ret <- list(data = dat, 
              data.preprocessed = dat.preprocessed, 
              formula = ifelse(params$adjust$perform, formula.covariates, ""))
  return(ret)
  
}

# Function to adjust p-values
adjust.pvals <- function(df, method) {
  
  df$p.adjusted <- p.adjust(df$p.value, method = method)
  df <- df %>%
    dplyr::arrange(., exposure, omic, p.adjusted)
  
  return(df)
  
}

# Function to plot p-values by exposure
plot.pvals <- function(df, threshold) {
  
  plt <- df %>% dplyr::group_by(omic) %>%
    ggplot2::ggplot(., mapping = aes(x = omic, 
                                     y = -log10(p.adjusted), 
                                     label = paste(exposure, "\n", 
                                                   round(p.adjusted, 3)))) +
    ggplot2::geom_point(color = dplyr::case_when(-log10(df[["p.adjusted"]]) >= threshold ~ "#d95f02",
                                                 TRUE ~ "#7570b3"),
                        alpha = 0.7) +
    ggplot2::ylab("-log10(p-values)") +
    ggplot2::theme(axis.text.x = element_text(angle = 90,
                                              vjust = 0.5, hjust = 0.5))
  
  # Add labels (p-values) to points only if significant results
  if (dim(subset(df, -log10(p.adjusted) >= threshold))[1] > 0) {
    df.subset <- subset(df, -log10(p.adjusted) >= threshold)
    plt <- plt + ggrepel::geom_label_repel(data = df.subset)
  }
  
  return(plt)
  
}

# Function to plot results of bootstrapping
plot.boots <- function(df, x = "exposure", y = "omic", path.save) {
  
  # Add `sign` to indicate whether an exposure is "significant" or not
  df$sign <- ifelse(sign(df$`.lower` * df$`.upper`) > 0, "s", "ns")
  df$exposure <- sapply(df$exposure, function(x) {strsplit(as.character(x), "_")[[1]][2]})
  
  # Plot the coefficient estimates as points, add bars corresponding to 
  # confidence intervals computed from the bootstraps.
  # Add labels only to "significant" exposures
  plt <- df %>%
    ggplot2::ggplot(., mapping = aes_string(x = x, y = ".estimate")) +
    ggplot2::geom_point(aes(color = sign), data = df) +
    ggplot2::geom_errorbar(aes(ymax = `.upper`, ymin = `.lower`)) +
    suppressWarnings(ggrepel::geom_label_repel(data = subset(df, sign == "s"), 
                                               aes_string(x, ".estimate", label = x), 
                                               position = position_jitter(width = 0.1, 
                                                                          height = 0.1), 
                                               max.overlaps = 50)) +
    ggplot2::facet_wrap(as.formula(paste("~", y)), ncol = 2)
  ggsave(filename = "plot_boots.png", device = "png", path = path.save, 
         width = 20, height = 70, units = "in", limitsize = FALSE)
  
  return(plt)
}

# Pipeline to fit models in case of no Bootstrapping
pipeline.single <- function(data, params, model, formula, omic.var) {
  mod <- list()
  
  # Fit model
  if (model == "glm") {
    mod$mod <- glm(formula = formula, family = params$family.glm, 
                   data = data)
    
  } else if (model == "regularized") {
    # glmnet does not support the use of a formula -> create dataframes 
    # based on provided formula
    Y <- data %>%
      dplyr::select(., strsplit(formula, "~")[[1]][1]) %>%
      as.matrix(.)
    X <- data %>%
      dplyr::select(., strsplit(strsplit(formula, "~")[[1]][2], "[+]")[[1]]) %>%
      # Convert characters and factors to dummy variables
      fastDummies::dummy_cols(.data = ., 
                              select_columns = NULL, 
                              remove_first_dummy = TRUE, remove_selected_columns = TRUE) %>%
      as.matrix(.)
    
    # Confounders not regularized
    penalty.factor <- rep(1, dim(X)[2])
    vars.to.keep <- which(!(colnames(X) %in% colnames(X[, grepl("log.hcp", 
                                                                colnames(X))])))
    penalty.factor[vars.to.keep] <- 0
    
    # Fit model
    mod$mod <- glmnet::cv.glmnet(x = X, y = Y, 
                                 nfolds = 10, parallel = TRUE, 
                                 family = params$regularization$family, 
                                 alpha = params$regularization$type, 
                                 standardize = FALSE, standardize.response = FALSE, 
                                 intercept = params$regularization$intercept, 
                                 penalty.factor = penalty.factor)
  } else if (model == "pls") {
    # Change default options for current session
    pls::pls.options(parallel = 4)
    
    # Fit model
    mod$mod <- pls::mvr(formula = as.formula(formula), ncomp = params$pls$ncomp, 
                        data = data, method = params$pls$method, 
                        scale = FALSE, center = FALSE, 
                        validation = "CV", jackknife = TRUE)
    
    return(mod)
  } # End if check type model
  
  # Tidy results
  if (params$adjust$perform == FALSE) {
    mod$mod <- broom::tidy(mod$mod)[2, ]
  } else {
    mod$mod <- broom::tidy(mod$mod) %>%
      dplyr::filter(stringr::str_starts(term, "log.hcp"))
  }
  mod$mod$omic <- omic.var
  
  return(mod)
}

# Helper function to fit one model to one sample of bootstrap
fit.model.split <- function(split, model, formula, params) {
  if (model == "glm") {
    fit <- glm(formula = formula, family = params$family.glm, 
               data = rsample::analysis(split))
  } else if (model == "regularized") {
    
  }
  
  return(fit)
}

# Pipeline to fit models in case of Bootstrapping
pipeline.bootstrapping <- function(data, params, model, formula, omic.var) {
  # Create samples
  boots <- rsample::bootstraps(data = data, 
                               times = params$boot$times, 
                               apparent = FALSE)
  
  # Fit one model for each sample
  boot_models <- boots %>%
    dplyr::mutate(model = purrr::map(splits, fit.model.split, 
                                     model = model, 
                                     formula = formula, 
                                     params = params), 
                  coef_info = purrr::map(model, broom::tidy))
  
  # Compute CI using percentile method for the exposures
  for (boot in 1:nrow(boot_models)) {
    tmp <- boot_models[boot, ]$coef_info[[1]] %>%
      dplyr::filter(stringr::str_starts(term, "log.hcp"))
    boot_models[boot, ]$coef_info <- list(tmp)
  }
  percentile_intervals <- suppressWarnings(rsample::int_pctl(boot_models, coef_info))
  
  # Tidy results
  percentile_intervals$omic <- omic.var
  
  return(percentile_intervals)
}

# Simple regression: 1 -omics ~ 1 exposure (+ covariates)
simple.oneToOne <- function(data, params, formula.covariates) {
  
  # Dataframe for all the results
  res <- data.frame()
  
  # Iterate over the -omics
  for (omic.var in colnames(data$omics)) {
    
    # Iterate over the exposures
    for (exposure in colnames(data$exposures)) {
      formula.mod <- ifelse(!params$adjust$perform | params$adjust$method == "residual", 
                            paste0(omic.var, "~", exposure), 
                            paste0(omic.var, "~", exposure, "+", formula.covariates))
      if (params$adjust$perform == FALSE | params$adjust$method == "residual") {
        data.mod <- cbind(data$omics, data$exposures)
      } else {
        data.mod <- cbind(data$omics, data$exposures, data$metadata)
      }
      
      # Check whether to perform bootstrapping on each model
      if (params$boot$perform == TRUE) {
        percentile_intervals <- pipeline.bootstrapping(data = data.mod, params = params, 
                                                       model = params$method.regression, 
                                                       formula = formula.mod, 
                                                       omic.var = omic.var)
        res <- rbind(res, percentile_intervals)
        
      } else {
        # In this case fit the model to the dataset only once
        mod <- pipeline.single(data = data.mod, params = params, 
                               model = params$method.regression, 
                               formula = formula.mod, 
                               omic.var = omic.var)
        res <- rbind(res, mod$mod)
      }
      
    } # End loop exposures
    
  } # End loop -omics
  
  res <- dplyr::rename(res, exposure = term)
  return(res)
  
}

# Simple regression: 1 -omics ~ all exposures (+ covariates)
simple.oneToMany <- function(data, params, formula.covariates) {
  
  # Dataframe for all the results
  res <- data.frame()
  # List containing fitted models PLS for all response variables
  mods.pls <- list()
  
  # Iterate over the -omics
  for (omic.var in colnames(data$omics)) {
    exposures <- paste(colnames(data$exposures), sep = "", collapse = "+")
    
    formula.mod <- ifelse(!params$adjust$perform | params$adjust$method == "residual", 
                          paste0(omic.var, "~", exposures), 
                          paste0(omic.var, "~", exposures, "+", formula.covariates))
    if (params$adjust$perform == FALSE | params$adjust$method == "residual") {
      data.mod <- cbind(data$omics, data$exposures)
    } else {
      data.mod <- cbind(data$omics, data$exposures, data$metadata)
    }
    
    # Check whether to perform bootstrapping on each model
    if (params$boot$perform == TRUE) {
      percentile_intervals <- pipeline.bootstrapping(data = data.mod, params = params, 
                                                     model = params$method.regression, 
                                                     formula = formula.mod, 
                                                     omic.var = omic.var)
      res <- rbind(res, percentile_intervals)
      
    } else {
      # In this case fit the model to the dataset only once
      #formula.mod <- "cbind(Sucrose, Proline.betaine) ~ log.hcp_detp_cadj + log.hcp_dmtp_cadj + log.hcp_dep_cadj"
      mod <- pipeline.single(data = data.mod, params = params, 
                             model = params$method.regression, 
                             formula = formula.mod, 
                             omic.var = omic.var)
      
      # Store fitted models PLS in order to print results
      if (params$method.regression == "pls") {
        mods.pls[[omic.var]] <- mod$mod
      } else {
        res <- rbind(res, mod$mod)
      }
    } # End if bootstrapping
    
  } # End loop -omics
  
  # Plot results
  if (params$method.regression == "pls") {
    plot.results.pls(mod = mods.pls, data = data.mod, 
                     formula = formula.mod, type = "single")
    return(mods.pls)
  } else {
    # In this case (right now (20/08/2021) only GLM models), returns dataframe with results
    res <- dplyr::rename(res, exposure = term)
    return(res)
  }
  
}

# Plot results fitted model from package PLS
plot.results.pls <- function(mod, data, formula, type) {
  # Argument `mod` in this case is a named list: one PLS model for each 
  # response variable
  
  ##### HELPER FUNCTIONS #####
  
  # Helper function to plot coefficients for 1 component of PLS
  plot.coefs.comp <- function(mod, data, formula, comp) {
    # Threshold to add error bars and labels (for predictors)
    threshold = 0.1
    
    # Compute SE
    se.all.models <- lapply(mod, 
                            function(x) {sqrt(pls::var.jack(x, ncomp = comp))}) %>%
      as.data.frame() %>% tibble::rownames_to_column(., var = "exposure") %>%
      tidyr::pivot_longer(!exposure, names_to = "omic", values_to = "se")
    
    mod <- mod %>%
      lapply(., function(x) { coef(x, comps = comp) }) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., var = "exposure") %>%
      tidyr::pivot_longer(!exposure, names_to = "omic", values_to = "coef")
    
    dat <- cbind(mod, se.all.models %>% dplyr::select(!c(exposure, omic)))
    dat$sign <- sign((dat$coef + dat$se) * (dat$coef - dat$se))
    
    plt <- ggplot2::ggplot(dat, mapping = aes(x = omic, y = coef, 
                                              label = exposure)) +
      ggplot2::geom_point(data = dat, mapping = aes(color = exposure)) +
      ggplot2::geom_errorbar(mapping = aes(ymin = coef - se, 
                                           ymax = coef + se, color = exposure), 
                             data = subset(dat, abs(coef) >= threshold & sign > 0)) +
      ggplot2::xlab("omic") + ggplot2::ylab("coefficient") +
      suppressWarnings(ggrepel::geom_label_repel(data = subset(dat, 
                                                               abs(coef) >= threshold & sign > 0), 
                                                 max.overlaps = 100, 
                                                 position = position_jitter())) +
      ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5, hjust = 0.5)) +
      ggplot2::ggtitle(paste0("Coefficients for component ", comp, "."))
    
    return(plt)
  }
  
  # Helper function to plot correlation circles
  plot.corr.circle <- function(mods) {
    radii <- c(sqrt(1/2), 1)
    df.segments <- data.frame(x1 = c(-1, 0), y1 = c(0, -1), 
                              x2 = c(1, 0), y2 = c(0, 1))
    
    for (mod.name in names(mods)) {
      mod <- mods[[mod.name]]
      
      scores.pls <- pls::scores(mod)[, 1:2] %>% as.data.frame()
      scores.x <- cor(model.matrix(mod), scores.pls) %>% as.data.frame()
      scores.x$label <- colnames(model.matrix(mod))
      scores.y <- cor(model.response(model.frame(mod)), scores.pls) %>%
        as.data.frame()
      scores.y <- tibble::rownames_to_column(scores.y, "label")
      scores.bind <- rbind(scores.x, scores.y) %>% as.data.frame()
      colnames(scores.bind) <- c("comp1", "comp2", "label")
      scores.bind$shape <- c(rep("predictor", dim(model.matrix(mod))[2]), 
                             rep("response", dim(scores.y)[1]))
      
      # Temporarily change variables' names to make them shorter
      labels.short <- scores.bind %>%
        dplyr::rowwise() %>%
        dplyr::transmute(label = ifelse(grepl("log.hcp", label), 
                                        strsplit(label, split = "_")[[1]][2], 
                                        label))
      scores.bind$label <- NULL
      scores.bind <- dplyr::bind_cols(scores.bind, label = labels.short)
      
      plt <- scores.bind %>%
        ggplot2::ggplot(data = ., mapping = aes(x = comp1, y = comp2, 
                                                label = label)) +
        ggplot2::geom_point(data = scores.bind, 
                            mapping = aes(color = label, 
                                          shape = shape, size = 2)) +
        ggplot2::coord_fixed(ratio = 1) +
        ggplot2::xlim(c(-1, 1)) + ggplot2::ylim(c(-1, 1)) +
        ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 1), n = 200, 
                             inherit.aes = FALSE, alpha = 0.1, size = 0.1) +
        ggplot2::geom_segment(mapping = aes(x = 0, y = 1, 
                                            xend = 0, yend = -1), 
                              alpha = 0.1, size = 0.1) +
        ggplot2::geom_segment(mapping = aes(x = -1, y = 0, 
                                            xend = 1, yend = 0), 
                              alpha = 0.1, size = 0.1) +
        ggplot2::geom_segment(data = subset(scores.bind, 
                                            abs(sqrt(comp1^2 + comp2^2)) >= 0.5), 
                              mapping = aes(x = 0, y = 0, 
                                            xend = comp1, yend = comp2, 
                                            color = label), 
                              arrow = arrow(length = unit(1/2, "picas"))) +
        suppressWarnings(ggrepel::geom_label_repel(data = subset(scores.bind, 
                                                                 abs(sqrt(comp1^2 + comp2^2)) >= 0.5), 
                                                   max.overlaps = 100)) +
        ggplot2::xlab("component 1") + ggplot2::ylab("component 2") +
        ggplot2::ggtitle(paste0("Correlation (>= 0.5) plot for ", mod.name, ".")) +
        ggplot2::theme(legend.position = "none")
      
      print(plt)
    } # End loop PLS models
    
  }
  
  # Helper function to plot X variance explained by each component
  plot.x.expl.var <- function(mods) {
    num.comp <- dim(pls::scores(mods[[1]]))[2]
    df.all <- data.frame()
    
    # Iterate over all PLS models and append results to dataframe
    for (mod.name in names(mods)) {
      mod <- mods[[mod.name]] %>% pls::explvar(.) %>% t() %>%
        as.data.frame() %>%
        `colnames<-`(c(paste0(rep("comp", num.comp), 1:num.comp))) %>%
        dplyr::mutate(., omic = mod.name)
      df.all <- rbind(df.all, mod)
    }
    
    # Plot
    plt <- df.all %>%
      tidyr::gather(., "comp", "var.expl", !omic) %>%
      ggplot2::ggplot(data = ., mapping = aes(x = reorder(omic, var.expl), 
                                              y = var.expl, 
                                              fill = factor(comp, 
                                                            levels = c(paste0(rep("comp", 
                                                                                  num.comp), 
                                                                              num.comp:1))))) +
      ggplot2::geom_bar(position = "stack", stat = "identity") +
      ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5, hjust = 0.5)) +
      ggplot2::xlab("omic") + ggplot2::ylab("percentage") +
      ggplot2::ggtitle("X variance explained by component.")
    
    return(plt)
    
  }
  
  # Helper function to plot X loadings by component
  plot.x.loadings <- function(mods) {
    num.comp <- dim(pls::scores(mods[[1]]))[2]
    
    for (mod.name in names(mod)) {
      mod <- mods[[mod.name]]$loadings[, 1:num.comp] %>% as.data.frame() %>%
        `colnames<-`(c(paste0(rep("comp", num.comp), 1:num.comp))) %>%
        tibble::rownames_to_column(., var = "exposure") %>%
        tidyr::pivot_longer(!exposure, names_to = "comp", values_to = "loading")
      
      plt <- mod %>%
        ggplot2::ggplot(data = ., mapping = aes(x = reorder(exposure, loading), 
                                                y = loading, color = comp)) +
        ggplot2::geom_line(data = mod, 
                           mapping = aes(group = comp)) +
        ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                  vjust = 0.5, hjust = 0.5)) +
        ggplot2::xlab("exposure") + ggplot2::ylab("loading") +
        ggplot2::ggtitle(paste0("Loadings of X by component for ", mod.name, "."))
      
      print(plt)
    } # End loop PLS models
    
  }
  
  # Helper function to plot X loadings by -omic variables
  plot.x.loadings.by.omic <- function(mods, variate) {
    df.all <- data.frame()
    
    # Iterate over all PLS models and append results to dataframe
    for (mod.name in names(mods)) {
      mod <- mods[[mod.name]]$loadings[, variate] %>% as.data.frame() %>%
        `colnames<-`(c(paste0("comp", variate))) %>%
        tibble::rownames_to_column(., var = "exposure") %>%
        tidyr::pivot_longer(!exposure, names_to = "comp", values_to = "loading") %>%
        dplyr::mutate(., omic = mod.name)
      
      df.all <- rbind(df.all, mod)
    }
    
    # Temporarily change variables' names to make them shorter
    labels.short <- df.all %>%
      dplyr::rowwise() %>%
      dplyr::transmute(label = ifelse(grepl("log.hcp", exposure), 
                                      strsplit(exposure, split = "_")[[1]][2], 
                                      exposure))
    df.all <- dplyr::bind_cols(df.all, label = labels.short)
    
    # Plot
    plt <- df.all %>%
      ggplot2::ggplot(data = ., mapping = aes(x = reorder(omic, loading), 
                                              y = loading, 
                                              color = exposure, 
                                              label = label)) +
      ggplot2::geom_point() +
      ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5, hjust = 0.5)) +
      suppressWarnings(ggrepel::geom_label_repel(data = subset(df.all, abs(loading) >= 5.0), 
                                                 max.overlaps = 100)) +
      ggplot2::xlab("omic") + ggplot2::ylab("loading") +
      ggplot2::ggtitle(paste0("Loadings by -omic of exposures for variate: ", variate, "."))
    
    print(plt)
    
  }
  
  # Helper function to plot Y loadings
  plot.y.loadings <- function(mods, variate) {
    df.all <- data.frame()
    
    # Iterate over all PLS models and append results to dataframe
    for (mod.name in names(mods)) {
      mod <- mods[[mod.name]]$Yloadings[, variate] %>% as.data.frame() %>%
        `colnames<-`(c(paste0("comp", variate))) %>%
        tibble::rownames_to_column(., var = "omic") %>%
        tidyr::pivot_longer(!omic, names_to = "comp", values_to = "loading") %>%
        dplyr::mutate(., omic = mod.name)
      
      df.all <- rbind(df.all, mod)
    }
    
    # Plot
    plt <- df.all %>%
      ggplot2::ggplot(data = ., mapping = aes(x = reorder(omic, loading), 
                                              y = loading, 
                                              color = comp, 
                                              label = omic)) +
      ggplot2::geom_point() +
      ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5, hjust = 0.5)) +
      # suppressWarnings(ggrepel::geom_label_repel(data = subset(df.all, abs(loading) >= 5.0), 
      #                                            max.overlaps = 100)) +
      ggplot2::xlab("omic") + ggplot2::ylab("loading") +
      ggplot2::ggtitle(paste0("Loadings of -omic variables for variate: ", variate, "."))
    
    print(plt)
    
  }
  
  # Helper function to plot scores
  plot.scores <- function(mods, data, var, col.metadata) {
    color.by <- data$metadata[, c(col.metadata)]
    
    # Iterate over fitted models and plot samples by variate 1 vs. variate 2
    for (mod.name in names(mods)) {
      if (var == "X") {
        mod <- mods[[mod.name]]$scores[, 1:2] %>% as.data.frame()
      } else if (var == "Y") {
        mod <- mods[[mod.name]]$Yscores[, 1:2] %>% as.data.frame()
      }
      mod <- cbind(mod, color.by)
      colnames(mod) <- c("comp1", "comp2", col.metadata)
      
      plt <- mod %>%
        ggplot2::ggplot(., mapping = aes(x = comp1, y = comp2)) +
        ggplot2::geom_point(mapping = aes(color = cohort.x)) +
        ggplot2::xlab("variate 1") + ggplot2::ylab("variate 2") +
        ggplot2::ggtitle(paste0("Variate 1 vs. variate 2 for samples (", var, 
                                ") - ", mod.name))
      
      print(plt)
    }
  }
  
  ##### PLOTTING #####
  
  # Plot coefficients
  # for (comp in 1:length(pls::compnames(mod[[1]]))) {
  #   # Each figure will contain the coefficients for all the exposures (on y-axis)
  #   # and all -omics (on x-axis), for one component of PLS
  #   plt.coefs <- plot.coefs.comp(mod, data, formula, comp)
  #   print(plt.coefs)
  # 
  # } # End loop components to plot coefficients
  
  # Plot of loadings for predictors
  # plot.x.loadings(mod)
  # plot.x.loadings.by.omic(mod, variate = 1)
  
  # Plot of loadings for responses
  # plot.y.loadings(mod, variate = 1)
  
  # Plot of X and Y scores (for samples)
  # plot.scores(mods = mod, data = data, 
  #             var = "X", col.metadata = "cohort.x")
  # plot.scores(mods = mod, data = data, 
  #             var = "Y", col.metadata = "cohort.x")
  
  # Correlation plot (circle)
  # plt.corr.circle <- plot.corr.circle(mod)
  
  # X variance explained by the components
  # plt.x.expl.var <- plot.x.expl.var(mod)
  # print(plt.x.expl.var)
  
}

##########
# TODO
# print(pls::selectNcomp(mod))
# print(pls::validationplot(mod, "R2, RMSEP, MSEP"))
# print(summary(mod))
##########

# Multi-Outcome regression: all -omics ~ 1 exposure (+ covariates)
multi.manyToOne <- function() {
  
}

# Multi-Outcome regression: all -omics ~ all exposures (+ covariates)
multi.manyToMany <- function() {
  
}
