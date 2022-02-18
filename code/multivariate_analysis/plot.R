# Author: Lorenzo Fabbri
# Script containing functions to plot results mixOmics

library(mixOmics)
library(tibble)
library(knitr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(labelled)
library(tidyr)
library(ggrepel)
library(gridExtra)
library(psych)

# Function to plot boxplots by categorical variable
boxplot.vars.by.group <- function(dat, group) {
  # Ungroup data to have one column with all variables
  dat <- tidyr::gather(dat, key = "var", value = "val", -group)
  
  gg.plot <- ggplot2::ggplot(dat, aes(x = var, y = val, color = group)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 1))
  print(gg.plot)
}

# Function to plot predicted vs. actual values
plot.predicted <- function(model) {
  X <- tibble::as_tibble(model$X)
  
  plots.tmp <- list()
  # Loop over the number of components of the model
  for (component in c(1)) {
  #for (component in 1:model$ncomp) {
    # Actual values
    features <- colnames(model$Y)
    Y <- model$Y %>% t() %>% as.data.frame() %>%
      dplyr::mutate(median.actual = apply(., MARGIN = 1, FUN = median), 
                    dev.actual = apply(., MARGIN = 1, FUN = mad)) %>%
      `rownames<-`(features) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::select(c(feature, median.actual, dev.actual)) %>%
      dplyr::arrange(feature)
    
    # Predicted values
    Y.hat <- predict(model, X)
    features <- colnames(Y.hat$predict[, , component])
    Y.hat$predict <- Y.hat$predict[, , component] %>%
      t() %>% as.data.frame() %>%
      dplyr::mutate(median.predicted = apply(., MARGIN = 1, FUN = median), 
                    dev.predicted = apply(., MARGIN = 1, FUN = mad)) %>%
      `rownames<-`(features) %>%
      tibble::rownames_to_column("feature") %>%
      dplyr::select(c(feature, median.predicted, dev.predicted)) %>%
      dplyr::arrange(feature)
    
    # Plot
    pred.plot <- ggplot2::ggplot(mapping = aes(x = Y$median.actual, 
                                               y = Y.hat$predict$median.predicted)) +
      ggplot2::geom_point(aes(color = colors()[component+1])) +
      ggplot2::geom_errorbar(data = Y.hat$predict,
                             aes(ymin = median.predicted - dev.predicted,
                                 ymax = median.predicted + dev.predicted), 
                             size = 0.1) +
      ggplot2::geom_errorbar(data = Y,
                             aes(xmin = median.actual - dev.actual,
                                 xmax = median.actual + dev.actual), 
                             size = 0.1) +
      ggplot2::geom_abline(slope = 1, intercept = 0) +
      # ggplot2::xlim(min(Y.hat$predict$median.predicted),
      #               max(Y.hat$predict$median.predicted)) +
      ggplot2::ggtitle(paste0("Predicted vs. Actual Values\nVariate: ", component)) +
      ggplot2::xlab("actual") + ggplot2::ylab("predicted") +
      ggplot2::theme(legend.position = "none")
    
    # Append plot to list
    plots.tmp <- append(plots.tmp, list(pred.plot))
  }
  
  pred.plots <- gridExtra::grid.arrange(grobs = plots.tmp, 
                                        #ncol = model$ncomp)
                                        ncol = 1)
  return(pred.plots)
}

# Function to plot PC_i vs. PC_(i+1) for subjects
plot.individuals <- function(model, metadata, color.by, is.adjusted, 
                             pc_i, pc_j) {
  variates <- model$variates
  variates.X <- as.data.frame(variates$X)
  variates.Y <- as.data.frame(variates$Y)
  
  # Convert indexes to strings to access variates
  pc_i <- paste0("comp", pc_i)
  pc_j <- paste0("comp", pc_j)
  
  if (!is.numeric(metadata[[color.by]])) {
    metadata$color.by <- as.factor(metadata$color.by)
    shape <- metadata[[color.by]]
    size <- 0.5
  } else {
    shape <- "1"
    size <- metadata[[color.by]]
  }
  
  if (is.adjusted$perform == FALSE | is.adjusted$method == "residual") {
    sex <- "e3_sex.x"
  } else { sex <- "e3_sex.x_female" }
  
  # X
  plot.indiv.X <- ggplot2::ggplot(data = variates.X, 
                                  mapping = aes_string(x = pc_i, y = pc_j)) +
    ggplot2::geom_point(aes(color = as.factor(metadata[[sex]]), 
                            shape = shape, size = size), 
                        alpha = 0.5) +
    ggplot2::ggtitle("Plot of Individuals for X")
  
  # Y
  plot.indiv.Y <- ggplot2::ggplot(data = variates.Y, 
                                  mapping = aes_string(x = pc_i, y = pc_j)) +
    ggplot2::geom_point(aes(color = as.factor(metadata[[sex]]), 
                            shape = shape, size = size), 
                        alpha = 0.5) +
    ggplot2::ggtitle("Plot of Individuals for Y")
  
  return(list(plot.indiv.X, plot.indiv.Y))
}

# Function to plot loadings of both X and Y per variate
plot.loadings <- function(model, model.code, variate, pls.params) {
  loadings.all <- model$loadings
  threshold.loading <- 0.1
  is.serum.metab <- strsplit(model.code, "_")[[1]][2] %>%
    stringr::str_sub(., 5, 5) == "2"
  title.y <- ifelse(is.serum.metab, 
                    paste0("Loadings for Y (>=", threshold.loading, ")"), 
                    "Loadings for Y")
  
  # If Serum Metabolome, filter results since too many variables
  if (is.serum.metab) {
    loadings.all$Y <- loadings.all$Y %>%
      data.frame() %>%
      dplyr::filter(abs(get(variate)) >= threshold.loading)
  }
  
  # X
  loadings.X <- loadings.all$X %>%
    as.data.frame() %>% rownames_to_column(var = "feature") %>%
    dplyr::select(c(feature, variate)) %>%
    ggplot2::ggplot(aes(x = reorder(feature, get(variate)), y = get(variate), 
                        fill = get(variate) > 0)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "feature", y = "loading", 
                  title = "Loadings for X", 
                  subtitle = paste0("Variate: ", variate, ".", " Model: ", model.code, sep = "")) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() + ggplot2::guides(fill = "none")
  
  # Y
  loadings.Y <- loadings.all$Y %>%
    as.data.frame() %>% rownames_to_column(var = "feature") %>%
    dplyr::select(c(feature, variate)) %>%
    ggplot2::ggplot(aes(x = reorder(feature, get(variate)), y = get(variate), 
                        fill = get(variate) > 0)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "feature", y = "loading", 
                  title = title.y, 
                  subtitle = paste0("Variate: ", variate, ".", " Model: ", model.code, sep = "")) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() + ggplot2::guides(fill = "none")
  
  plots <- gridExtra::grid.arrange(grobs = list(loadings.X, loadings.Y), ncol = 2)
  
  if (pls.params$plotting$save == TRUE) {
    ggplot2::ggsave(filename = paste0("loadings_", variate, "_", model.code, ".png", sep = ""), 
                    plot = plots, 
                    device = "png", 
                    path = paste0(pls.params$path.save, "plots_by_time/"), 
                    dpi = 320)
  }
  
}

# Helper function to make correlation circles (`plotVar` in mixOmics)
plot.cor.circles <- function(model, model.code, c_i, c_j) {
  # This function uses Pearson correlation coefficients
  
  radii <- c(sqrt(1/2), 1)
  df.segments <- data.frame(x1 = c(-1, 0), y1 = c(0, -1), 
                            x2 = c(1, 0), y2 = c(0, 1))
  
  if (class(model) == "mixo_pls") {
    cordX.1 <- cor(model$X, model$variates$X[, c(c_i, c_j)], use = "pairwise")
    cordX.2 <- cor(model$Y, model$variates$X[, c(c_i, c_j)], use = "pairwise")
    
    # Test of significance correlation
    corr.test.1 <- psych::corr.test(model$X, model$variates$X[, c(c_i, c_j)])
    corr.test.2 <- psych::corr.test(model$Y, model$variates$X[, c(c_i, c_j)])
    
  } else if (class(model) == "mixo_spls") {
    ll.X <- unique(unlist(lapply(c(c_i, c_j), function(x) {
      mixOmics::selectVar(model, comp = x)$X$name
    })))
    ll.Y <- unique(unlist(lapply(c(c_i, c_j), function(x) {
      mixOmics::selectVar(model, comp = x)$Y$name
    })))
    
    cordX.1 <- cor(model$X[, colnames(model$X) %in% ll.X, drop = FALSE], 
                   model$variates$X[, c(c_i, c_j)], use = "pairwise")
    cordX.2 <- cor(model$Y[, colnames(model$Y) %in% ll.Y, drop = FALSE], 
                   model$variates$X[, c(c_i, c_j)], use = "pairwise")
  }
  
  scores.bind <- rbind(cordX.1, cordX.2) %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "label")
  scores.bind$shape <- c(rep("predictor", dim(cordX.1)[1]), 
                         rep("response", dim(cordX.2)[1]))
  
  # Temporarily change variables' names to make them shorter
  scores.bind <- scores.bind %>%
    dplyr::mutate(score = round(abs(sqrt(.[[2]]^2 + .[[3]]^2)), 2))
  
  labels.short <- scores.bind %>%
    dplyr::rowwise() %>%
    dplyr::transmute(label = ifelse(grepl("log.hcp", label), 
                                    paste0(strsplit(label, split = "_")[[1]][2], 
                                           "\n", score), 
                                    paste0(label, "\n", score)))
  scores.bind$label <- NULL
  scores.bind <- dplyr::bind_cols(scores.bind, label = labels.short)
  
  # Plot
  threshold <- 0.3
  comp1 <- paste0("comp", c_i)
  comp2 <- paste0("comp", c_j)
  
  plt <- scores.bind %>%
    ggplot2::ggplot(data = ., mapping = aes_string(x = paste0("comp", c_i), 
                                                   y = paste0("comp", c_j), 
                                                   label = "label")) +
    ggplot2::geom_point(data = scores.bind, 
                        mapping = aes(color = label, shape = shape,
                                      size = 2)) +
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
                                        abs(sqrt(scores.bind[[1]]^2 + scores.bind[[2]]^2)) >= threshold), 
                          mapping = aes_string(x = 0, y = 0, 
                                               xend = comp1, yend = comp2), 
                          arrow = arrow(length = unit(1/2, "picas"))) +
    suppressWarnings(ggrepel::geom_label_repel(data = subset(scores.bind, 
                                                             abs(sqrt(scores.bind[[1]]^2 + scores.bind[[2]]^2)) >= threshold), 
                                               max.overlaps = 100)) +
    ggplot2::xlab(comp1) + ggplot2::ylab(comp2) +
    ggplot2::ggtitle(paste0("Correlation circle with a threshold of ", threshold, 
                            ".\nModel: ", model.code)) +
    ggplot2::theme(legend.position = "none")
  
  return(plt)
  
}

# Function to plot internal methods of `mixOmics`
plot.internal <- function(model, perf.obj, metadata, 
                          pls.params) {
  
  model.code = pls.params$model.code
  color.by.indiv = pls.params$color.by.indiv
  is.adjusted = pls.params$perform.adj
  
  # Variable plot
  if (model$ncomp > 1) {
    plts.circles <- list()
    for (i in 1:(model$ncomp - 1)) {
      # mixOmics::plotVar(model, legend = TRUE, comp = c(i, i + 1), cutoff = 0.2, 
      #                   title = "Variable Plot", style = "ggplot2")
      
      plt.tmp <- plot.cor.circles(model = model, model.code = model.code, 
                                  c_i = i, c_j = i + 1)
      print(plt.tmp)
      if (pls.params$plotting$save == TRUE) {
        file.name <- paste0("corr_comp", i, "comp", i + 1, "_")
        ggplot2::ggsave(filename = paste0(file.name, model.code, ".png"), 
                        device = "png", 
                        path = paste0(pls.params$path.save, "plots_by_time/"), 
                        dpi = 320)
      }
      plts.circles[[i]] <- plt.tmp
    }
    #gridExtra::grid.arrange(grobs = plts.circles, ncol = 1)
  }
  
  # Loadings plot
  for (comp in 1:model$ncomp) {
    #mixOmics::plotLoadings(model, comp = comp)
    plot.loadings(model = model, model.code = model.code, 
                  variate = paste0("comp", comp), pls.params = pls.params)
  }
  
  # Individuals plot
  #mixOmics::plotIndiv(object = model)
  if (model$ncomp > 1) {
    for (i in 1:(model$ncomp - 1)) {
      plot.indiv <- plot.individuals(model = model, metadata = metadata, 
                                     color.by = color.by.indiv, 
                                     is.adjusted = is.adjusted, 
                                     pc_i = i, pc_j = i + 1)
      gridExtra::grid.arrange(grobs = plot.indiv, 
                              ncol = 2)
    }
  }
  
  # Clustered Image Map
  for (comp in 1:model$ncomp) {
    mixOmics::cim(model, comp = comp, margins = c(7, 7),
                  title = paste0("CIM for component: ", comp))
  }
  
}

# Function to plot Variate i vs. Variate (i+1) and color by <Q2> of i
plot.variates.Q2 <- function(model, labels, pls.params, 
                             metric.summary, dat, threshold.Q2, 
                             pc_i, pc_j) {
  model.code <- pls.params$model.code
  loadings <- model$loadings
  
  # if (pc_i > 1) {
  #   threshold.Q2 <- 0.05
  # }
  
  # Re-arrange summary metric
  data.tmp <- metric.summary %>%
    dplyr::select(-c(sd, min, max)) %>%
    spread(key = comp, value = mean) %>%
    dplyr::select(feature, 
                  as.character(pc_i), as.character(pc_j))
  colnames(data.tmp) <- c("feature", 
                          paste0("comp", pc_i, ".Q2"), 
                          paste0("comp", pc_j, ".Q2"))
  data.tmp$feature <- as.character(data.tmp$feature)
  
  # Add <Q2> to summary metric
  loadings.dat <- loadings[[dat]] %>%
    as.data.frame() %>%
    dplyr::select(paste0("comp", pc_i), paste0("comp", pc_j)) %>%
    tibble::rownames_to_column(var = "feature")
  
  if (dat == "Y") {
    old.nrow <- nrow(loadings.dat)
    loadings.dat <- loadings.dat %>%
      dplyr::select("feature", 
                    paste0("comp", pc_i), paste0("comp", pc_j)) %>%
      dplyr::full_join(., data.tmp)
    if (old.nrow != nrow(loadings.dat)) {stop(call. = TRUE)}
  } else {
    loadings.dat[[paste0("comp", pc_i, ".Q2")]] <- threshold.Q2
    loadings.dat[[paste0("comp", pc_j, ".Q2")]] <- threshold.Q2
  }
  
  # Add family exposures (i.e., labels for features)
  if (dat == "X") {
    old.nrow <- nrow(loadings.dat)
    loadings.dat <- loadings.dat %>%
      dplyr::inner_join(labels)
    if (old.nrow != nrow(loadings.dat)) {stop(call. = TRUE)}
    loadings.dat[[paste0("comp", pc_i, ".Q2")]] <- loadings.dat$family
  }
  
  # Plot
  title <- paste0("Plot of Variate ", pc_i, " vs. Variate ", pc_j, 
                  "\nColored by average Q2 of Variate ", pc_i, "\n")
  title <- paste0(title, dat)
  loadings.plot <- ggplot2::ggplot(loadings.dat, 
                                   aes_string(x = paste0("comp", pc_i), 
                                              y = paste0("comp", pc_j), 
                                              label = "feature")) +
    ggplot2::geom_point(aes_string(color = paste0("comp", pc_i, ".Q2"))) +
    ggrepel::geom_text_repel(data = loadings.dat %>% dplyr::mutate(feature = ifelse(abs(!!sym(paste0("comp", 
                                                                                                     ifelse(dat == "X", pc_j, pc_i), 
                                                                                                     ".Q2"))) >= threshold.Q2, 
                                                                                    feature, "")), 
                             size = 3.5, 
                             box.padding = 0.25, max.overlaps = 30) +
    ggplot2::labs(title = title, 
                  subtitle = paste0("Model: ", model.code))
  file.name <- paste0("comp", pc_i, "comp", pc_j, "_")
  
  if (pls.params$plotting$save == TRUE) {
    ggplot2::ggsave(filename = paste0(file.name, dat, "_", model.code, ".png"), 
                    device = "png", 
                    path = paste0(pls.params$path.save, "plots_by_time/"), 
                    dpi = 320)
  }
  
  return(loadings.plot)
}

# Helper function to plot a single metric from mixOmics object
plot.metric <- function(model, measures, var, pls.params, labels) {
  model.code <- pls.params$model.code
  
  old.var <- var
  if (var == "R2.comp") {
    var <- stringr::str_split(var, "[.]")[[1]][1]
  }
  metric              <- measures[[var]]
  metric.summary      <- metric$summary
  metric.summary$comp <- as.factor(metric.summary$comp)
  metric.values       <- metric$values
  metric.values$comp  <- as.factor(metric.values$comp)
  
  if (old.var == "R2.comp") {
    pc.plots.tmp <- list()
    for (pc in levels(metric.summary$comp)) {
      # Re-compute each time since I modify the same variable in each loop
      metric.summary      <- metric$summary
      metric.summary$comp <- as.factor(metric.summary$comp)
      
      # For this value we add more information to the plot
      values <- metric.values %>%
        dplyr::filter(comp == pc) %>%
        dplyr::group_by(feature) %>%
        dplyr::summarise(min = min(value), 
                         max = max(value))

      nrow.old <- nrow(metric.summary %>% dplyr::filter(comp == pc))
      metric.summary <- dplyr::inner_join(metric.summary %>% dplyr::filter(comp == pc), 
                                          values[, c("feature", "min", "max")])
      if (nrow(metric.summary) != nrow.old) {
        stop(call. = TRUE)
      }
      
      # If Serum Metabolome, filter results since too many variables
      threshold.r2.comp <- 0.03
      is.serum.metab <- strsplit(model.code, "_")[[1]][2] %>%
        stringr::str_sub(., 5, 5) == "2"
      if (is.serum.metab) {
        metric.summary <- metric.summary %>%
          data.frame() %>%
          dplyr::filter(mean >= threshold.r2.comp)
      } # End filtering Serum Metabolome
      
      pc.plot.tmp <- metric.summary %>%
        dplyr::arrange(comp, desc(mean)) %>%
        dplyr::mutate(feature = factor(feature, levels = feature)) %>%
        ggplot2::ggplot(aes(x = feature, y = mean, fill = comp)) +
        ggplot2::geom_col(position = position_dodge(0.6), width = 0.6) +
        ggplot2::theme(legend.position = "none", axis.title.y = element_blank()) +
        ggplot2::labs(title = paste0(var, " Summary: ", pc), 
                      subtitle = paste0("Model: ", model.code)) +
        ggplot2::coord_flip()
      
      if (pls.params$plotting$save == TRUE) {
        ggplot2::ggsave(filename = paste0("r2_summary_comp", pc, "_", model.code, ".png", sep = ""), 
                        device = "png", 
                        path = paste0(pls.params$path.save, "plots_by_time/"), 
                        dpi = 320)
      }
      
      # Print also table since difficult to add labels to plot
      # print(metric.summary %>%
      #         dplyr::mutate_at(vars(mean, min, max, sd), list(~ round(., 2))) %>%
      #         dplyr::arrange(desc(mean)))
      
      pc.plots.tmp <- append(pc.plots.tmp, list(pc.plot.tmp))
    }
    metric.summary.plot <- ggpubr::ggarrange(plotlist = pc.plots.tmp, ncol = 3)
  } else if (old.var == "Q2") {
    # For this value we add more information to the plot
    values <- metric$values %>%
      dplyr::group_by(feature, comp) %>%
      dplyr::summarise(min = min(value), 
                       max = max(value))
    values$comp <- as.factor(values$comp)
    
    nrow.old <- nrow(metric.summary)
    metric.summary <- dplyr::inner_join(metric.summary, 
                                        values[, c("feature", "comp", "min", "max")])
    if (nrow(metric.summary) != nrow.old) {
      stop(call. = TRUE)
    }
    
    metric.summary <- metric.summary %>%
      dplyr::arrange(comp, desc(mean)) %>%
      dplyr::mutate(feature = factor(feature, levels = unique(feature)))
    col.labs <- rep(FALSE, nrow(metric.summary))
    col.labs[metric.summary$mean > 0.0] <- TRUE
    
    # If Serum Metabolome, filter results since too many variables
    metric.summary.old <- metric.summary
    is.serum.metab <- strsplit(model.code, "_")[[1]][2] %>%
      stringr::str_sub(., 5, 5) == "2"
    if (is.serum.metab) {
      metric.summary <- metric.summary %>%
        data.frame() %>%
        dplyr::filter(mean >= -0.01)
    } # End filtering Serum Metabolome
    
    metric.summary.plot <- ggplot2::ggplot(metric.summary, 
                                           aes(x = feature, 
                                               y = mean, 
                                               fill = comp)) +
      ggplot2::geom_point(aes(colour = comp)) +
      ggplot2::labs(title = paste0(var, " Summary"), 
                    subtitle = paste0("Model: ", model.code)) +
      ggplot2::coord_flip()
    metric.summary.plot <- metric.summary.plot +
      ggplot2::geom_text(data = subset(metric.summary, mean > 0.0), 
                         mapping = aes(x = feature, 
                                       y = mean, 
                                       label = round(mean, 3)), 
                         position = position_jitter(width  = 0.07, 
                                                    height = 0.05), 
                         size = 3) +
      ggplot2::geom_vline(data = metric.summary[col.labs, ], 
                          aes(xintercept = feature), 
                          linetype = "dashed", color = "blue", size = 0.2) +
      ggplot2::geom_hline(data = metric.summary, 
                          aes(yintercept = 0.0975), 
                          linetype = "dashed", color = "blue")
    
    if (pls.params$plotting$save == TRUE) {
      ggplot2::ggsave(filename = paste0("q2_summary_", model.code, ".png"), 
                      device = "png", 
                      path = paste0(pls.params$path.save, "plots_by_time/"), 
                      dpi = 320)
    }
    
    # Print also table since difficult to add labels to plot
    # print(metric.summary %>% dplyr::group_by(comp) %>%
    #         dplyr::mutate_at(vars(mean, min, max, sd), list(~ round(., 2))) %>%
    #         dplyr::arrange(desc(mean)))
    
    ## Plot of variate 1 vs. variate 2 for both X and Y
    ## Color by <Q2> of each exposure/omic
    
    if (model$ncomp > 1) {
      for (i in 1:(model$ncomp - 1)) {
        # Y (-omics)
        is.serum.metab <- strsplit(pls.params$model.code, "_")[[1]][2] %>%
          stringr::str_sub(., 5, 5) == "2"
        threshold.Q2.Y <- ifelse(is.serum.metab, 0.03, -1)
        loadings.Y.plot <- plot.variates.Q2(model = model, labels = labels, 
                                            pls.params = pls.params, 
                                            metric.summary = metric.summary.old, 
                                            dat = "Y", threshold.Q2 = threshold.Q2.Y, 
                                            pc_i = i, pc_j = i + 1)
        # X (exposures)
        loadings.X.plot <- plot.variates.Q2(model = model, labels = labels, 
                                            pls.params = pls.params, 
                                            metric.summary = metric.summary.old, 
                                            dat = "X", threshold.Q2 = -1, 
                                            pc_i = i, pc_j = i + 1)
        
        gridExtra::grid.arrange(loadings.X.plot, loadings.Y.plot, 
                                ncol = 2)
      }
    }
    
  } else if (old.var == "Q2.total") {
    # For this value we add more information to the plot
    values <- metric$values %>%
      dplyr::group_by(comp) %>%
      dplyr::summarise(min = min(value), 
                       max = max(value))
    values$comp <- as.factor(values$comp)
    
    nrow.old <- nrow(metric.summary)
    metric.summary <- dplyr::inner_join(metric.summary, 
                                        values[, c("comp", "min", "max")])
    if (nrow(metric.summary) != nrow.old) { stop(call. = TRUE) }
    
    # Write information to file
    readr::write_lines(metric.summary, 
                       file = paste0(pls.params$path.save, "tables/", 
                                     model.code, ".txt", collapse = ""), 
                       append = TRUE)
    
    metric.summary.plot <- ggplot2::ggplot(metric.summary, aes(x = comp, 
                                                               y = mean)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(data = metric.summary, 
                             aes(ymin = mean - sd, 
                                 ymax = mean + sd), 
                             width = 0.1) +
      ggplot2::geom_hline(yintercept = 0.0975, 
                          linetype = "dashed", color = "blue") +
      ggplot2::geom_text(mapping = aes(comp, mean, 
                                       label = paste0("(", round(mean, 3), ", ", 
                                                      round(min, 2), ", ", 
                                                      round(max, 2), ")")), 
                         position = position_jitter(height = 0.05)) +
      ggplot2::labs(title = paste0(var, " Summary"), 
                    subtitle = paste0("Model: ", model.code))
    
    if (pls.params$plotting$save == TRUE) {
      ggplot2::ggsave(filename = paste0("q2_total_", model.code, ".png"), 
                      device = "png", 
                      path = paste0(pls.params$path.save, "plots_by_time/"), 
                      dpi = 320)
    }
    
  } else {
    # If Serum Metabolome, filter results since too many variables
    threshold.r2 <- 0.03
    is.serum.metab <- strsplit(model.code, "_")[[1]][2] %>%
      stringr::str_sub(., 5, 5) == "2"
    if (is.serum.metab & old.var == "R2") {
      metric.summary <- metric.summary %>%
        data.frame() %>%
        dplyr::filter(mean >= threshold.r2)
    } # End filtering Serum Metabolome
    
    metric.summary.plot <- ggplot2::ggplot(metric.summary %>%
                                             dplyr::arrange(comp, desc(mean)) %>%
                                             dplyr::mutate(feature = factor(feature, levels = unique(feature))), 
                                           aes(x = feature, y = mean, fill = comp)) +
      ggplot2::geom_point(aes(colour = comp)) +
      ggplot2::ggtitle(paste0(var, " Summary")) +
      ggplot2::coord_flip()
  }
  
  return(metric.summary.plot)
}

# Function to visualize performance metrics
plot.perf <- function(model, perf.obj, data, pls.params) {
  model.code <- pls.params$model.code
  
  ## Model
  prop_expl_var <- tibble::as_tibble(model$prop_expl_var)
  print(prop_expl_var)
  
  # Write information to file
  readr::write_lines(prop_expl_var, 
                     file = paste0(pls.params$path.save, "tables/", 
                                   model.code, ".txt", collapse = ""), 
                     append = TRUE)
    
  ## Perf
  # Labels exposures (family)
  labels.exp <- data$exposures %>%
    labelled::var_label() %>%
    unlist() %>% as.data.frame() %>%
    tibble::rownames_to_column('feature') %>%
    dplyr::rename(., family = `.`)
  
  # Temporarily add family to exposures not present in results Proteome
  problematic_exposures <- list("log.hcp_detp_cadj" = "OP Pesticides", 
                                "log.hcp_mepa_cadj" = "Phenols")
  for (el in names(problematic_exposures)) {
    if (!(el %in% labels.exp$feature)) {
      labels.exp <- labels.exp %>%
        tibble::add_row(feature = el, family = problematic_exposures[[el]])
    }
  }
  
  measures <- perf.obj$measures
  vars <- c("MSEP", # MSEP (Mean Square Error Prediction for each Y variable)
            #"RMSEP", # RMSEP (Root Mean Square Error or Prediction for each Y variable)
            "R2", "R2.comp", # R2 (R2 for each Y variable)
            "Q2", # Q2 (Q2 for each Y variable)
            "Q2.total" # Q2 total (# iterations CV * # components)
            #"RSS", # RSS (Residual Sum-of-Squares for each Y variable)
            #"PRESS", # PRESS (Predicted Residual Error Sum-of-Squares for each Y variable)
            #"cor.tpred", "cor.upred", # cor.tpred & cor.upred (correlation between actual and predicted components
                                      # of X (t) and Y (u) -> # iterations CV * # components)
            #"RSS.tpred", "RSS.upred" # RSS.tpred & RSS.upred (Residual Sum-of-Squares between actual and predicted
                                     # components of X (t) and Y (u))
            )
  for (var in vars) {
    print(plot.metric(model = model, measures = measures,
                      var = var, pls.params = pls.params,
                      labels = labels.exp))
  }
  
  # VIP
  vip.pls <- mixOmics::vip(model)
  colnames(labels.exp) <- c("exposure", "family")
  
  plot.vip.by.exp <- vip.pls %>%
    as.data.frame() %>%
    tibble::rownames_to_column('exposure') %>%
    dplyr::full_join(labels.exp, by = "exposure") %>%
    tibble::as_tibble() %>%
    tidyr::gather(key = 'key', value = 'value', 
                  -exposure, -family) %>%
    ggplot2::ggplot(aes(x = key, y = value, fill = exposure, group = value)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap( ~ family) +
    ggplot2::labs(title = "VIP by exposure family", 
                  subtitle = paste0("Model: ", model.code)) +
    ggplot2::guides(fill = guide_legend(ncol = 2))
  print(plot.vip.by.exp)
        
  # Stability (only for sPLS)
  if ("stability.X" %in% names(perf.obj$features)) {
    stab.x <- perf.obj$features$stability.X
    stab.y <- perf.obj$features$stability.Y
    
    # Exposures
    stab.x.plot <- stab.x %>%
      as.data.frame() %>% tibble::rownames_to_column('exposure') %>%
      dplyr::full_join(labels.exp, by = "exposure") %>%
      tidyr::gather(key = 'key', value = 'value', 
                    -exposure, -family) %>%
      ggplot2::ggplot(aes(x = key, y = value, fill = exposure)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap( ~ family) +
      ggplot2::ggtitle("Stability for X") +
      ggplot2::guides(fill = guide_legend(ncol = 2))
    print(stab.x.plot)
    
    # Omics
    stab.y.plot <- stab.y %>%
      as.data.frame() %>% tibble::rownames_to_column('omics') %>%
      tidyr::gather(key = 'key', value = 'value', 
                    -omics) %>%
      ggplot2::ggplot(aes(x = key, y = value, fill = omics)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::ggtitle("Stability for Y") +
      ggplot2::guides(fill = guide_legend(ncol = 2))
    print(stab.y.plot)
    
    # Print table by component showing top contributors
    stab.y.table <- stab.y %>%
      as.data.frame() %>% tibble::rownames_to_column('omics') %>%
      tidyr::gather(key = 'key', value = 'value', 
                    -omics) %>%
      dplyr::group_split(key) %>% lapply(., arrange, desc(value))
    print(stab.y.table)
  }
  
}
