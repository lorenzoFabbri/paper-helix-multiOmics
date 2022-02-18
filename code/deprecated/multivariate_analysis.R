# Author: Lorenzo Fabbri
# Functions to perform Multivariate Analyses

###############################################################################
################################## DEPRECATED #################################
###############################################################################

library(mixOmics)
library(dplyr)
library(ggpubr)
library(tibble)
library(stringr)
library(labelled)
library(tidyr)
library(magrittr)
library(ggh4x)
library(scales)
library(tidytext)
library(parallel)

# Function to perform all analyses with simple (s)PLS
perform.simple <- function(omics,
                           exposures,
                           metadata,
                           pls.params,
                           is.permute = FALSE,
                           iter.seed = NULL) {
  # Unpack parameters for (s)PLS
  num.components <- pls.params$num.components
  num.vars       <- pls.params$num.vars
  is.sparse      <- pls.params$is.sparse
  include.bmi    <- pls.params$include.bmi
  perform.adj    <- pls.params$perform.adj
  scale.perf     <- pls.params$scale.perf
  log.transform  <- pls.params$log.transform
  
  # Run simple (s)PLS
  res.spls.simple <- run.spls.simple(
    omics,
    exposures,
    metadata,
    num.components,
    num.vars,
    is.sparse,
    include.bmi,
    perform.adj,
    scale.perf,
    log.transform,
    is.permute
  )
  perf.spls.simple <- res.spls.simple[[3]]
  tune.spls.simple <- res.spls.simple[[1]]
  res.spls.simple  <- res.spls.simple[[2]]
  
  if (is.permute == FALSE) {
    # Plot results internal to mixOmics
    plot.pls(res.spls.simple, perf.spls.simple)
    
    # Plot metrics extracted from results
    plot.pls.metrics(res.spls.simple, perf.spls.simple, exposures)
  }
  
  return(list(tune.spls.simple, res.spls.simple, perf.spls.simple))
}

# Helper function to handle missing values
handle.missing.values <- function(dat, type) {
  if (type == 'remove') {
    dat <- dat %>%
      tidyr::drop_na()
  }
  
  return(dat)
}

# Helper function to pre-process datasets
change.type <- function(dat, col, new.type) {
  dat[[col]] %<>% new.type
  return(dat)
}

preprocess.data <- function(exposome,
                            omic,
                            metadata,
                            scale.perf,
                            log.transform,
                            is.permute = FALSE) {
  # Drop rows containing missing values
  exposome <- handle.missing.values(exposome, 'remove')
  omic     <- handle.missing.values(omic, 'remove')
  
  exposome <- change.type(exposome, "SubjectID", as.character)
  omic     <- change.type(omic, "SubjectID", as.character)
  metadata <- change.type(metadata, "SubjectID", as.character)
  
  # Check that datasets contain same subjects
  not.present <- dplyr::anti_join(exposome, omic)
  if (dim(not.present)[1] > 0) {
    cat('Removing missing subjects in Exposome...\n')
    exposome <- exposome %>%
      dplyr::filter(!(SubjectID %in% not.present$SubjectID))
    metadata <- metadata %>%
      dplyr::filter(!(SubjectID %in% not.present$SubjectID))
  }
  not.present <- dplyr::anti_join(omic, exposome)
  if (dim(not.present)[1] > 0) {
    cat('Removing missing subjects in Omics...\n')
    omic <- omic %>%
      dplyr::filter(!(SubjectID %in% not.present$SubjectID))
    metadata <- metadata %>%
      dplyr::filter(!(SubjectID %in% not.present$SubjectID))
  }
  
  # Check order rows is same in all datasets
  cat('Checking order of rows...\n')
  if (is.permute == FALSE) {
    print(sum(exposome$SubjectID == omic$SubjectID) == dim(exposome)[1])
    print(sum(exposome$SubjectID == metadata$SubjectID) == dim(exposome)[1])
  } else {
    # Shuffle rows datasets
    exposome[] <- lapply(exposome, sample)
    omic[]     <- lapply(omic, sample)
    metadata[] <- lapply(metadata, sample)
    
    print(sum(exposome$SubjectID == omic$SubjectID) != dim(exposome)[1])
    print(sum(exposome$SubjectID == metadata$SubjectID) != dim(exposome)[1])
  }
  
  # log-transform variables if necessary
  if (log.transform) {
    # Check if "log" is present in column names
    if (length(omic %>% dplyr::select(dplyr::contains("log")))) {
      cat("Some columns seem to be already log-transformed.")
    }
    
    omic <- omic %>%
      dplyr::transmute(dplyr::across(where(is.numeric), log2,
                                     .names = "log.{.col}"))
  }
  
  # Scale and center data by cohort
  debug.scaling <- FALSE
  if (scale.perf) {
    # Check if data contain NaN
    if (sum(is.na(exposome)) != 0) {
      print('MISSING VALUES IN EXPOSOME!')
    }
    if (sum(is.na(omic)) != 0) {
      print('MISSING VALUES IN OMICS!')
    }
    
    scale.tibble <- function(x) {
      (x - mean(x)) / sd(x)
    }
    
    scale.by.cohort <- function(dat) {
      # Extract code cohorts
      cohort <- stringr::str_sub(dat$HelixID, start = 1, end = 3)
      dat <- tibble::add_column(dat, cohort)
      
      dat.normalized <- dat %>%
        dplyr::select(-contains('desc')) %>%
        group_by(cohort) %>%
        dplyr::mutate_at(vars(-c(HelixID, SubjectID, cohort)), scale.tibble)
      
      return(dat.normalized)
    }
    
    exposome <- scale.by.cohort(exposome)
    omic     <- scale.by.cohort(omic)
    
    if (debug.scaling) {
      boxplot.vars.byGroup <- function(dat, title) {
        dat <- tidyr::gather(dat, key = "key", value = "val",
                             3:(dim(dat)[2] - 1))
        
        gg.plot <-
          ggplot2::ggplot(dat, aes(
            x = key,
            y = val,
            color = cohort
          )) +
          geom_boxplot() +
          theme(axis.text.x = element_text(
            angle = 60,
            vjust = 0.5,
            hjust = 1
          ))
        print(gg.plot)
        # ggplot2::ggsave(filename = paste0("/home/lorenzo/Downloads/", title, ".png"),
        #                 width = 20, height = 20)
      }
      
      # Check distribution variables after scaling
      boxplot.vars.byGroup(exposome, 'exp')
      boxplot.vars.byGroup(omic, 'omic')
    }
  }
  
  # Drop columns containing metadata and categorical variables (e.g., *_cdesc)
  omic <- omic %>%
    dplyr::select(-c(HelixID, SubjectID)) %>%
    dplyr::select(where(is.numeric)) %>%
    ungroup()
  exposome <- exposome %>%
    dplyr::select(-c(HelixID, SubjectID)) %>%
    dplyr::select(where(is.numeric)) %>%
    ungroup()
  
  if (scale.perf) {
    omic <- omic %>%
      dplyr::select(-c(cohort)) %>%
      ungroup()
    exposome <- exposome %>%
      dplyr::select(-c(cohort)) %>%
      ungroup()
  }
  
  # Check density -omics after scaling by cohort
  if (scale.perf) {
    if (debug.scaling) {
      tmp.data <- omic %>%
        tidyr::gather() %>%
        dplyr::filter(key != "cohort") %>%
        dplyr::ungroup() %>%
        ggplot(aes(value)) +
        geom_density() +
        stat_theodensity(color = "blue") +
        facet_wrap( ~ key, scales = "free") +
        ggtitle(paste('Distribution variables after cohort-scaling'))
      path.plot <-
        paste0("/home/lorenzo/Downloads/",
               "_distAfterCohortScaling",
               ".png")
      ggsave(path.plot,
             plot = tmp.data,
             width = 20,
             height = 20)
    }
  }
  
  return(list(exposome, omic, metadata))
}

# Function to perform 2-stage residual-outcome regression
# (since cannot adjust for covariates in xPLS, yet...)
perform.2sr <- function(outcome, covariates) {
  res.model <- lm(as.matrix(outcome) ~ as.matrix(covariates))
  predicted <- predict(res.model, covariates, type = 'response')
  
  # Deflation step
  deflated <- outcome - predicted
  
  return(deflated)
}

# Simple (s)PLS: omic~exposures
run.spls.simple <- function(omic,
                            exposome,
                            metadata,
                            num.components,
                            num.vars,
                            is.sparse,
                            include.bmi,
                            perform.adj,
                            scale.perf,
                            log.transform,
                            is.permute) {
  # Input:
  # num.components refers to the latent factors of PLS
  # num.vars       refers to the number of variables to keep in sPLS
  
  # Pre-process datasets
  tmp <- preprocess.data(
    exposome,
    omic,
    metadata,
    scale.perf = scale.perf,
    log.transform = log.transform,
    is.permute
  )
  exposome <- tmp[[1]]
  omic     <- tmp[[2]]
  metadata <- tmp[[3]]
  
  # Write datasets to file for further debugging, if necessary
  readr::write_csv(exposome, file = "./data/exposome.csv",
                   col_names = TRUE)
  readr::write_csv(omic, file = "./data/omics.csv",
                   col_names = TRUE)
  readr::write_csv(metadata, file = "./data/metadata.csv",
                   col_names = TRUE)
  
  # Adjust for covariates using GLM
  covariates <- colnames(metadata)
  if (perform.adj) {
    covariates <- c(
      'cohort',
      'ethn_dummy2',
      'age_sample_years',
      'e3_sex',
      '',
      '',
      '',
      '',
      'Urine.SamplingType'
    )
    
    # TEMPORARY! Just for testing only zBMI
    covariates <- c('zBMI')
    print("IMPORTANT: ONLY zBMI IS USED NOW!")
    if (!include.bmi) {
      covariates <- c()
    }
    
    omic <- perform.2sr(omic, metadata[, covariates])
  }
  
  # Perform (s)PLS (by default sparsity on both X *and* Y)
  run.model <- function(exposome,
                        omic,
                        is.sparse,
                        is.test,
                        ncomp,
                        keepX,
                        keepY,
                        max.iter) {
    if (is.sparse == FALSE) {
      spls.results <- mixOmics::pls(
        X = exposome,
        Y = omic,
        ncomp = ncomp,
        mode = 'regression',
        scale = FALSE,
        max.iter = max.iter
      )
    } else {
      if (is.test == TRUE) {
        # If we do not care about sparsity
        spls.results <- mixOmics::spls(
          X = exposome,
          Y = omic,
          ncomp = ncomp,
          mode = 'regression',
          scale = FALSE,
          max.iter = max.iter
        )
      } else {
        spls.results <- mixOmics::spls(
          X = exposome,
          Y = omic,
          ncomp = ncomp,
          keepX = keepX,
          keepY = keepY,
          mode = 'regression',
          scale = FALSE,
          max.iter = max.iter
        )
      }
    }
    
    return(spls.results)
  }
  
  pipeline <- function(exposome,
                       omic,
                       is.sparse,
                       num.components) {
    # Parameters and hyper-parameters
    cv.nrepeat <- 5
    cv.folds   <- 5
    max.iter   <- 500
    
    # Run model (model)
    spls.results <- run.model(
      exposome = exposome,
      omic = omic,
      is.sparse = is.sparse,
      is.test = TRUE,
      ncomp = num.components,
      keepX = NULL,
      keepY = NULL,
      max.iter = max.iter
    )
    
    # Check performances (perf) to select ncomp
    spls.perf <- mixOmics::perf(
      spls.results,
      validation = 'Mfold',
      folds = cv.folds,
      nrepeat = cv.nrepeat,
      criterion = 'all',
      progressBar = FALSE
    )
    print(plot(spls.perf, criterion = 'Q2.total'))
    
    # Tune hyper-parameters (tune) for keepX and/or keepY
    if (is.sparse == TRUE) {
      # Apparently "pls" is not accepted in "tune"...
      spls.tune <-
        mixOmics::tune(
          method = "spls",
          mode = "regression",
          X = exposome,
          Y = omic,
          ncomp = num.components,
          test.keepY = c(1, 3, 5, 7, ncol(exposome)),
          test.keepX = c(1, 3, 5, 7, ncol(omic)),
          nrepeat = cv.nrepeat,
          validation = "Mfold",
          folds = cv.folds,
          measure = "cor",
          scale = FALSE,
          center = FALSE,
          near.zero.var = FALSE,
          cpus = parallel::detectCores(),
          max.iter = max.iter
        )
      print(plot(spls.tune, measure = 'cor'))
      
      # Check performances (perf) with chosen (hyper-)parameters
      spls.results <- run.model(
        exposome = exposome,
        omic = omic,
        is.sparse = is.sparse,
        is.test = FALSE,
        ncomp = spls.tune$choice.ncomp,
        keepX = spls.tune$choice.keepX,
        keepY = spls.tune$choice.keepY,
        max.iter = max.iter
      )
      spls.perf <- mixOmics::perf(
        spls.results,
        validation = 'Mfold',
        folds = cv.folds,
        nrepeat = cv.nrepeat,
        criterion = 'all',
        progressBar = FALSE
      )
      
      cat("\n==================================================\n")
      cat("Dimension of Exposome (X):\n")
      print(dim(exposome))
      cat("Dimension of Omics (Y):\n")
      print(dim(omic))
      cat("Summary of (s)PLS:\n")
      cat("Number of components selected:\n")
      print(spls.tune$choice.ncomp)
      cat("Number of keepX selected:\n")
      print(spls.tune$choice.keepX)
      cat("Number of keepY selected:\n")
      print(spls.tune$choice.keepY)
      cat("\n==================================================\n")
      cat("\n")
    } else {
      # PLS
      spls.tune <- NULL
      
      cat("\n==================================================\n")
      cat("Dimension of Exposome (X):\n")
      print(dim(exposome))
      cat("Dimension of Omics (Y):\n")
      print(dim(omic))
      cat("\n==================================================\n")
      cat("\n")
    }
    
    return(list(spls.tune, spls.results, spls.perf))
  }
  
  res.pipeline <-
    pipeline(exposome, omic, is.sparse, num.components)
  spls.tune    <- res.pipeline[[1]]
  spls.results <- res.pipeline[[2]]
  spls.perf    <- res.pipeline[[3]]
  
  return(list(spls.tune, spls.results, spls.perf))
}

# Function to plot results from PLS in mixOmics
plot.pls <- function(pls.results, spls.perf) {
  # Variable plot
  if (pls.results$ncomp > 1) {
    mixOmics::plotVar(pls.results, legend = TRUE,
                      #cutoff = 0.5,
                      title = 'Variable plot')
  }
  
  # Loadings plot
  for (comp in pls.results$ncomp) {
    mixOmics::plotLoadings(pls.results, comp = comp)
  }
  
  # Clustered Image Map to show correlation structure between blocks
  # marg <- 7
  # for (comp in pls.results$ncomp) {
  #   mixOmics::cim(
  #     pls.results,
  #     comp = comp,
  #     margins = c(marg, marg),
  #     title = paste('CIM for component: ', comp)
  #   )
  # }
}

# Function to print selected results of (s)PLS model
plot.pls.metrics <- function(mod, perf, exposome) {
  ## MODEL
  # prop_expl_var (proportion of variance explained by each variate
  #               divided by the total variance in the data)
  prop_expl_var <- as_tibble(mod$prop_expl_var)
  knitr::kable(prop_expl_var, caption = "Proportion of variance explained")
  print(prop_expl_var)
  
  ## PERF (measures -> values & summary)
  measures <- perf$measures
  
  plot.metric.internal <- function(measures, var) {
    old.var <- var
    if (var == "R2.comp") {
      var <- str_split(var, "[.]")[[1]][1]
    }
    metric <- measures[[var]]
    metric.summary <- metric$summary
    metric.summary$comp <- as.factor(metric.summary$comp)
    
    if (old.var == "R2.comp") {
      pc.plots.tmp <- list()
      for (pc in levels(metric.summary$comp)) {
        pc.plot.tmp <- metric.summary %>%
          filter(comp == pc) %>%
          ggplot(aes(
            x = reorder(feature,-mean),
            y = mean,
            fill = comp
          )) +
          geom_col(position = position_dodge(0.6), width = 0.6) +
          theme(legend.position = "none",
                axis.title.y = element_blank()) +
          ggtitle(paste(var, " Summary: ", pc)) +
          coord_flip()
        
        pc.plots.tmp <- append(pc.plots.tmp,
                               list(pc.plot.tmp))
      }
      
      metric.summary.plot <-
        ggpubr::ggarrange(plotlist = pc.plots.tmp, ncol = 3)
      
    } else if (old.var == "Q2") {
      metric.summary <- arrange(metric.summary,-mean)
      col.labs <- rep(FALSE, nrow(metric.summary))
      col.labs[metric.summary$mean > 0.0] <- TRUE
      
      metric.summary.plot <- ggplot(metric.summary,
                                    aes(
                                      x = reorder(feature,-mean),
                                      y = mean,
                                      fill = comp
                                    )) +
        geom_point(aes(colour = comp)) +
        ggtitle(paste(var, " Summary")) +
        coord_flip()
      metric.summary.plot <- metric.summary.plot +
        geom_text(
          data = subset(metric.summary, mean > 0.0),
          mapping = aes(
            x = reorder(feature,-mean),
            y = mean,
            label = round(mean, 3)
          ),
          position = position_jitter(width = 0.07,
                                     height = 0.05),
          size = 3
        ) +
        #theme(axis.text.y = element_text(color = col.labs)) +
        geom_vline(
          data = metric.summary[col.labs,],
          aes(xintercept = feature),
          # Rule of thumb
          linetype = "dashed",
          color = "blue",
          size = 0.2
        ) +
        geom_hline(
          data = metric.summary,
          aes(yintercept = 0.0975),
          # Rule of thumb
          linetype = "dashed",
          color = "blue"
        )
    } else if (old.var == "Q2.total") {
      metric.summary.plot <- ggplot(metric.summary,
                                    aes(x = comp, y = mean)) +
        geom_point() +
        geom_errorbar(data = metric.summary,
                      aes(ymin = mean - sd, ymax = mean + sd),
                      width = 0.1) +
        geom_hline(yintercept = 0.0975,
                   # Rule of thumb
                   linetype = "dashed",
                   color = "blue") +
        geom_text(
          mapping = aes(comp, mean,
                        label = round(mean, 3)),
          position = position_jitter(height = 0.05)
        ) +
        ggtitle(paste(var, " Summary"))
    } else {
      metric.summary.plot <- ggplot(metric.summary,
                                    aes(
                                      x = reorder(feature,-mean),
                                      y = mean,
                                      fill = comp
                                    )) +
        geom_point(aes(colour = comp)) +
        ggtitle(paste(var, " Summary")) +
        coord_flip()
    }
    
    return(metric.summary.plot)
  }
  
  # MSEP (Mean Square Error Prediction for each Y variable)
  print(plot.metric.internal(measures = measures, var = "MSEP"))
  
  # RMSEP (Root Mean Square Error or Prediction for each Y variable)
  print(plot.metric.internal(measures = measures, var = "RMSEP"))
  
  # R2 (R2 for each Y variable)
  print(plot.metric.internal(measures = measures, var = "R2"))
  print(plot.metric.internal(measures = measures, var = "R2.comp"))
  
  # Q2 (Q2 for each Y variable)
  print(plot.metric.internal(measures = measures, var = "Q2"))
  
  # Q2 total (# iterations CV * # components)
  print(plot.metric.internal(measures = measures, var = "Q2.total"))
  
  # RSS (Residual Sum-of-Squares for each Y variable)
  print(plot.metric.internal(measures = measures, var = "RSS"))
  
  # PRESS (Predicted Residual Error Sum-of-Squares for each Y variable)
  print(plot.metric.internal(measures = measures, var = "PRESS"))
  
  # cor.tpred & cor.upred (correlation between actual and predicted components
  #                       of X (t) and Y (u) -> # iterations CV * # components)
  cor.t <-
    plot.metric.internal(measures = measures, var = "cor.tpred")
  cor.u <-
    plot.metric.internal(measures = measures, var = "cor.upred")
  print(ggarrange(cor.t, cor.u))
  
  # RSS.tpred & RSS.upred (Residual Sum-of-Squares between actual and predicted
  #                       components of X (t) and Y (u))
  rss.t <-
    plot.metric.internal(measures = measures, var = "RSS.tpred")
  rss.u <-
    plot.metric.internal(measures = measures, var = "RSS.upred")
  print(ggarrange(rss.t, rss.u))
  
  # VIP
  vip.pls <- mixOmics::vip(mod)
  print(plot(vip.pls))
  
  labels.exp <- exposome %>%  # Labels exposures (family)
    var_label() %>%
    unlist() %>% as.data.frame() %>%
    rownames_to_column('exposure') %>% rename(., family = `.`)
  
  plot.vip.by.exp <- vip.pls %>%
    as.data.frame() %>%
    rownames_to_column('exposure') %>%
    full_join(labels.exp, by = "exposure") %>%
    as_tibble() %>%
    gather(key = 'key', value = 'value',-exposure,-family) %>%
    ggplot(aes(x = key, y = value, fill = exposure)) +
    geom_col(position = "dodge") +
    facet_wrap( ~ family) +
    ggtitle("VIP by exposure") +
    guides(fill = guide_legend(ncol = 2))
  print(plot.vip.by.exp)
  
  # Stability (only for sPLS)
  if ("stability.X" %in% names(perf$features)) {
    stab.x <- perf$features$stability.X
    stab.y <- perf$features$stability.Y
    
    stab.x.plot <- stab.x %>%  # Exposures
      as.data.frame() %>% rownames_to_column('exposure') %>%
      full_join(labels.exp, by = "exposure") %>%
      gather(key = 'key', value = 'value',-exposure,-family) %>%
      ggplot(aes(x = key, y = value, fill = exposure)) +
      geom_col(position = "dodge") +
      facet_wrap( ~ family) +
      ggtitle("Stability for X") +
      guides(fill = guide_legend(ncol = 2))
    print(stab.x.plot)
    
    stab.y.plot <- stab.y %>%  # Omics
      as.data.frame() %>% rownames_to_column('omics') %>%
      gather(key = 'key', value = 'value',-omics) %>%
      ggplot(aes(x = key, y = value, fill = omics)) +
      geom_col(position = "dodge") +
      ggtitle("Stability for Y") +
      guides(fill = guide_legend(ncol = 2))
    print(stab.y.plot)
    # Print also table by component showing top contributors
    stab.y.table <- stab.y %>%
      as.data.frame() %>% rownames_to_column('omics') %>%
      gather(key = 'key', value = 'value',-omics) %>%
      group_split(key) %>% lapply(., arrange, desc(value)) %>%
      print()
  }
}

# Function to perform validation by permutation
validate.permutation <- function(omic,
                                 exposome,
                                 metadata,
                                 pls.params,
                                 num.iters) {
  # Perform permutation num.iters times
  for (iter in 1:num.iters) {
    cat(paste("\nIteration number: ", iter, "\n"))
    res <- perform.simple(
      omics = omic,
      exposures = exposome,
      metadata = metadata,
      pls.params = pls.params,
      is.permute = TRUE,
      iter.seed = iter
    )
    tune <- res[[1]]
    res  <- res[[2]]
    perf <- res[[3]]
  }
  
}

# Multivariate sPLS: omics~exposures

# Multivariate & Multilevel sPLS: ?~?
