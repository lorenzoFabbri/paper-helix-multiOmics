# Script for generating images for paper

library(tidyverse)
library(tidygraph)

source("code/multivariate_analysis/dictionaries.R")

##### Function for QQ plot of p-values from results EWAS
qq.pval.ewas <- function(path) {
  
  ewas.t1 <- list.files(path = path, pattern = "^t1")
  ewas.t2 <- list.files(path = path, pattern = "^t2")
  ewas.t1.df <- list()
  ewas.t2.df <- list()
  for (idx in seq_along(ewas.t1)) {
    df1 <- read.csv(file = paste0(path, ewas.t1[[idx]])) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(fdr = p.adjust(p.value, "fdr")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(cpg = strsplit(cpg, split = "_me")[[1]][1])
    df2 <- read.csv(file = paste0(path, ewas.t2[[idx]])) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(fdr = p.adjust(p.value, "fdr")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(cpg = strsplit(cpg, split = "_me")[[1]][1])
    mol <- ewas.t1[[idx]] %>% strsplit(., split = "\\.") %>% .[[1]] %>% .[2] %>%
      strsplit(., split = "_") %>% .[[1]] %>% .[1]
    ewas.t1.df[[mol]] <- df1
    ewas.t2.df[[mol]] <- df2
  } # End loop over files results EWAS
  
  # Helper function for QQ plot of p-values
  # Largely copied from https://github.com/cran/QCEWAS/blob/master/R/script_v12-2_package.R
  qq.plot <- function(pvals, mol, lambda) {
    qq.expected <- sort(-log10(ppoints(length(pvals))))
    qq.observed <- sort(-log10(pvals))
    
    qq.exp.min <- qq.expected[1]
    qq.exp.max <- qq.expected[length(qq.expected)]
    qq.obs.min <- qq.observed[1]
    qq.obs.max <- qq.observed[length(qq.observed)]
    
    temp <- (1:length(pvals))
    i1000 <- c(1, (1:1000) * floor(length(pvals) / 1000), length(pvals))
    qq.band.upper <- sort(-log10(qbeta(1 - 0.05 / 2, temp, length(pvals) - temp + 1)))[i1000]
    qq.band.lower <- sort(-log10(qbeta(0.05 / 2, temp, length(pvals) - temp + 1)))[i1000]
    qq.expected.1000 <- qq.expected[i1000]
    
    ggplot2::ggplot() +
      ggplot2::geom_point(mapping = aes(x = qq.expected, y = qq.observed)) +
      ggplot2::geom_abline(colour = "red", show.legend = FALSE) +
      ggplot2::coord_cartesian(xlim = c(0, qq.exp.max), 
                               ylim = c(0, qq.obs.max)) +
      ggplot2::geom_polygon(mapping = aes(
        x = c(qq.expected.1000, rev(qq.expected.1000)), 
        y = c(qq.band.upper, rev(qq.expected.1000)), 
        alpha = 0.1
      ), show.legend = FALSE) +
      ggplot2::geom_polygon(mapping = aes(
        x = c(qq.expected.1000, rev(qq.expected.1000)), 
        y = c(qq.band.lower, rev(qq.expected.1000)), 
        alpha = 0.1
      ), show.legend = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Expected -log10(p-value)", 
                    y = "Observed -log10(p-value)", 
                    title = toupper(mol), 
                    subtitle = paste0("Lambda: ", round(lambda, 3)))
  }
  
  # T1
  list.qqplots <- list()
  for (mol in names(ewas.t1.df)) {
    lambda <- QCEWAS::P_lambda(ewas.t1.df[[mol]]$p.value)
    plt <- qq.plot(pvals = ewas.t1.df[[mol]]$p.value, 
                   mol = mol, lambda = lambda)
    list.qqplots <- append(list.qqplots, list(plt))
  }
  a.grob <- gridExtra::arrangeGrob(grobs = list.qqplots, ncol = 2)
  ggplot2::ggsave(filename = "qqplots_ewasT1.png", path = "results/images/", 
                  height = 12, width = 12, plot = a.grob)
  
  # T2
  list.qqplots <- list()
  for (mol in names(ewas.t2.df)) {
    lambda <- QCEWAS::P_lambda(ewas.t2.df[[mol]]$p.value)
    plt <- qq.plot(pvals = ewas.t2.df[[mol]]$p.value, 
                   mol = mol, lambda = lambda)
    list.qqplots <- append(list.qqplots, list(plt))
  }
  a.grob <- gridExtra::arrangeGrob(grobs = list.qqplots, ncol = 2)
  ggplot2::ggsave(filename = "qqplots_ewasT2.png", path = "results/images/", 
                  height = 12, width = 12, plot = a.grob)
  
}

##### Function to visualize common significant CpG sites by chemical
visualize.common.cpgs <- function(path, threshold) {
  ewas.t1 <- list.files(path = path, pattern = "^t1")
  ewas.t2 <- list.files(path = path, pattern = "^t2")
  
  ewas.t1.df <- list()
  ewas.t2.df <- list()
  for (idx in seq_along(ewas.t1)) {
    df1 <- read.csv(file = paste0(path, ewas.t1[[idx]])) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(fdr = p.adjust(p.value, "fdr"))
    df2 <- read.csv(file = paste0(path, ewas.t2[[idx]])) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(fdr = p.adjust(p.value, "fdr"))
    mol <- ewas.t1[[idx]] %>% strsplit(., split = "\\.") %>% .[[1]] %>% .[2] %>%
      strsplit(., split = "_") %>% .[[1]] %>% .[1]
    ewas.t1.df[[mol]] <- df1
    ewas.t2.df[[mol]] <- df2
  }
  
  list.venn.plts <- list()
  for (mol in names(ewas.t1.df)) {
    df1.signif <- ewas.t1.df[[mol]] %>%
      dplyr::filter(fdr <= threshold)
    df2.signif <- ewas.t2.df[[mol]] %>%
      dplyr::filter(fdr <= threshold)
    
    if (nrow(df1.signif) > 0 & nrow(df2.signif) > 0) {
      common.sites <- base::intersect(df1.signif$cpg, df2.signif$cpg)
      dat <- list(
        t1 = df1.signif$cpg, 
        t2 = df2.signif$cpg
      )
      
      venn <- ggVennDiagram::Venn(dat) %>% ggVennDiagram::process_data()
      
      venn.plt <- ggplot2::ggplot() +
        ggplot2::geom_sf(aes(fill = count), data = venn@region, 
                         show.legend = FALSE) +
        ggplot2::geom_sf(aes(color = id), size = 1, color = "white", 
                         data = venn@setEdge, 
                         show.legend = FALSE) +
        ggplot2::geom_sf_text(aes(label = name), data = venn@setLabel) +
        ggplot2::geom_sf_label(aes(label = count), 
                               family = "serif", size = 5, alpha = 0.5, 
                               data = venn@region) +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(toupper(mol), "\n", 
                                     "Threshold FDR: ", threshold))
      list.venn.plts[[mol]] <- venn.plt
    } else {
      common.sites <- "none"
    }
  } # End loop over molecules
  
  a.grob <- gridExtra::arrangeGrob(grobs = list.venn.plts, ncol = 2)
  ggplot2::ggsave(filename = "venn_molsCpG.png", path = "results/images/", 
                  height = 12, width = 12, plot = a.grob)
}

##### Function to load results bootstrapping merged networks and visualize intervals
plot.bootstrapping.nets <- function(path) {
  merged.nets <- list.files(path = path, 
                            pattern = "merged_net")
  
  net.all <- list()
  df.all <- list()
  for (i in seq_along(merged.nets)) {
    file.name <- merged.nets[[i]]
    net <- unname(get(load(paste0(path, file.name)))) %>% .[[1]] %>% .$net
    net.all[[i]] <- net
    
    node.attributes <- net %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
    edge.attributes <- net %>% tidygraph::activate(., what = "edges") %>%
      tidygraph::as_tibble()
    dim.old <- dim(edge.attributes)[1]
    df <- dplyr::inner_join(node.attributes, edge.attributes, 
                            by = c("idx" = "from"))
    df <- dplyr::inner_join(node.attributes, df, 
                            by = c("idx" = "to"))
    df <- df %>%
      dplyr::select(-c(group.x, idx, group.y, idx.y, type.x, type.y))
    if (dim(df)[1] != dim(edge.attributes)[1]) { stop(call. = TRUE) }
    df.all[[i]] <- df
  }
  
  # Find merged network of bootstraps of merged networks
  df.1 <- purrr::reduce(df.all, dplyr::inner_join, 
                        by = c("name.x", "name.y", 
                               "label.x", "label.y", 
                               "direction"))
  df.2 <- purrr::reduce(df.all, dplyr::inner_join, 
                        by = c("name.x" = "name.y", "name.y" = "name.x", 
                               "label.x" = "label.y", "label.y" = "label.x", 
                               "direction"))
  df <- dplyr::bind_rows(df.1, df.2)
  
  metrics <- c("pcor.x", "pval.x", "qval.x", "pcor.y", "pval.y", "qval.y")
  res.summaries <- list()
  for (met in metrics) {
    df.tmp <- df %>%
      dplyr::select(c(name.x, name.y, label.x, label.y, 
                      tidyselect::starts_with(met))) %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "var") %>%
      tibble::as_tibble()
    summ <- psych::describe(df.tmp[5:nrow(df.tmp), 2:ncol(df.tmp)] %>%
                              dplyr::mutate_all(as.numeric)) %>%
      data.frame() %>%
      dplyr::select(c(mean, sd, median, min, max, range)) %>%
      t() %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "var") %>%
      dplyr::mutate(var = paste0(var, "_", met))
    res.summaries <- append(res.summaries, 
                            list(rbind(df.tmp[1:4, ], summ)))
  }
  # Tibble where each column is an edge and each row (from the 5th) is
  # as numeric summary (e.g., mean, sd, ...)
  res.summaries <- dplyr::bind_rows(res.summaries) %>%
    dplyr::distinct()
  colnames(res.summaries) <- c("variable", paste0("E", 
                                                  1:(ncol(res.summaries)-1)))
  
  # Plot distribution of metrics of interest
  metrics <- c("pcor.x", "pcor.y")
  all.plots <- list(pcor.x = list(), pcor.y = list())
  for (metric in metrics) {
    df.tmp <- df %>%
      dplyr::select(c(name.x, name.y, 
                      tidyselect::starts_with(metric)))
    for (idx in 1:nrow(df.tmp)) {
      lab <- paste0(df.tmp[idx, 1], "-", df.tmp[idx, 2])
      plt <- df.tmp[idx, 3:ncol(df.tmp)] %>%
        c() %>% unname() %>% unlist() %>%
        as.data.frame() %>%
        ggplot2::ggplot() +
        ggplot2::geom_density(mapping = aes(x = `.`)) +
        ggplot2::geom_vline(mapping = aes(xintercept = median(`.`), 
                                          color = "median"), 
                            show.legend = TRUE) +
        ggplot2::labs(title = paste0(lab, "\n", 
                                     metric)) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_manual(name = "statistics", 
                                    values = c(median = "red")) +
        ggplot2::theme(axis.title.x = element_blank(), 
                       plot.title = element_text(size = 10))
      all.plots[[metric]] <- append(all.plots[[metric]], 
                                    list(plt))
    }
  }
  path.save <- "results/images/"
  grid.t1 <- gridExtra::grid.arrange(grobs = all.plots$pcor.x %>% unname(), ncol = 2, 
                                     nrow = length(all.plots$pcor.x) / 2)
  ggplot2::ggsave(filename = paste0(path.save, "pcor_t1_boot.pdf"), 
                  plot = grid.t1, 
                  dpi = 720/2, 
                  width = 7, height = 15)
  grid.t2 <- gridExtra::grid.arrange(grobs = all.plots$pcor.y %>% unname(), ncol = 2, 
                                     nrow = length(all.plots$pcor.x) / 2)
  ggplot2::ggsave(filename = paste0(path.save, "pcor_t2_boot.pdf"), 
                  plot = grid.t2, 
                  dpi = 720/2, 
                  width = 7, height = 15)
  
  # Compute statistics from bootstrapped merged networks
  df.props <- list()
  for (i in seq_along(net.all)) {
    net <- net.all[[i]]
    num.edges <- igraph::ecount(net)
    num.nodes <- igraph::vcount(net)
    num.components <- igraph::components(net)$no
    net.density <- igraph::graph.density(net)
    
    df.props[[i]] <- tibble::tibble(
      iteration = i, 
      num.edges = num.edges, num.nodes = num.nodes, 
      num.components = num.components, 
      net.density = net.density
      )
  }
  df.props.table <- reduce(df.props, dplyr::bind_rows) %>%
    gtsummary::tbl_summary() %>%
    gtsummary::bold_labels() %>%
    gtsummary::as_gt() %>% gt::gtsave(., 
                                      filename = "results/images/boot_desc.png", 
                                      zoom = 4, expand = 7)
}

##### Function to plot change in networks' metrics by time point
plot.change.time.metrics <- function(path.res, metrics) {
  load("../data/intermediate_res_ggm/processed_ggms.RData")
  
  net.1 <- processed.ggms[[1]]$net
  net.2 <- processed.ggms[[2]]$net
  
  tmp <- dplyr::inner_join(net.1 %>% as_tibble(), 
                           net.2 %>% as_tibble(), 
                           by = "name") %>%
    dplyr::mutate(change = (.[[metrics[2]]] - .[[metrics[1]]]) / .[[metrics[1]]]) %>%
    dplyr::arrange(desc(change)) %>%
    dplyr::top_n(n = 100) %>%
    dplyr::select(c(name, metrics[1], metrics[2])) %>%
    tidyr::gather(metric, value, -name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(name = strsplit(name, "_")[[1]][1])
  
  plt <- ggpubr::ggbarplot(
    data = tmp, x = "name", y = "value", fill = "metric", 
    x.text.angle = -45, 
    position = ggplot2::position_dodge(), 
    ggtheme = ggpubr::theme_pubclean()
  )
  plot(plt)
}

# Helper function to plot correlation matrix for Omics
internal.plotCorr <- function(data, data.feats, save.key) {
  path.save <- "results/images/"
  
  colnames(data) <- colnames(data) %>%
    lapply(., function(x) {
      strsplit(x, "_") %>% .[[1]] %>% .[1]
    })
  name.variables <- ifelse(grepl("prot", save.key), "Prot_ID", "Rvar")
  data <- data[data.feats[[name.variables]]]
  
  heatmap.omic <- data %>%
    cor() %>%
    ggcorrplot::ggcorrplot(type = "upper", show.diag = TRUE, 
                           legend.title = "cor", 
                           tl.cex = 4, tl.srt = 90)
  ggplot2::ggsave(paste0(path.save, save.key, ".png"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
}

##### Function to plot heatmaps of chemicals and -omics w/ additional info
plot.heatmaps <- function(time.point) {
  path.meta <- "data/"
  path.data <- "../data/"
  path.save <- "results/images/"
  
  # Exposures
  ## By chemical class but in same figure
  exposome <- readRDS(file = paste0(path.data, "exposome_1", 
                                    ifelse(time.point == 1, "A", "B"))) %>%
    .$data %>% tibble::as_tibble() %>%
    dplyr::select(-c(SampleID))
  colnames(exposome) <- colnames(exposome) %>%
    lapply(., function(x) {
      strsplit(x, "log.") %>% .[[1]] %>% .[2] %>%
        strsplit(., "_") %>% .[[1]] %>% .[1] %>%
        toupper(.)
    })
  all.exposures <- dict.exposure.groups() %>%
    stack() %>% tibble::as_tibble() %>%
    `colnames<-`(c("exposure", "class")) %>%
    dplyr::mutate(class = as.character(class)) %>%
    dplyr::mutate(class = substr(class, 1, nchar(class) - 1)) %>%
    dplyr::mutate(exposure = toupper(exposure))
  exposome <- exposome[all.exposures$exposure]
  
  cols.exp <- map.char.to.aes()[2][[1]] %>% as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble() %>%
    `colnames<-`(c("class", "color"))
  all.exposures <- all.exposures %>%
    dplyr::inner_join(cols.exp)
  
  heatmap.exp <- exposome %>% scale() %>%
    cor() %>%
    ggcorrplot::ggcorrplot(lab = TRUE, type = "upper", show.diag = TRUE, 
                           legend.title = "cor") +
    ggplot2::theme(axis.text.x = element_text(colour = all.exposures$color), 
                   axis.text.y = element_text(colour = all.exposures$color))
  ggplot2::ggsave(paste0(path.save, "corrPlot_exps", time.point, ".png"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
  
  # Omics
  ## Serum Metabolome
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                  ifelse(time.point == 1, 
                                         "A", "B")))$data %>%
    dplyr::select(tidyselect::ends_with("_ms")) %>% scale() %>%
    tibble::as_tibble()
  metabs.feats <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                        ifelse(time.point == 1, 
                                               "A", "B")))$feature.data
  internal.plotCorr(metabs, metabs.feats, 
                    paste0("corrPlot_metabs", time.point))
  
  ## Urinary Metabolome
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                  ifelse(time.point == 1, 
                                         "A", "B")))$data %>%
    dplyr::select(tidyselect::ends_with("_mu")) %>% scale() %>%
    tibble::as_tibble()
  metabu.feats <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                        ifelse(time.point == 1, 
                                               "A", "B")))$feature.data
  internal.plotCorr(metabu, metabu.feats, 
                    paste0("corrPlot_metabu", time.point))
  
  ## Proteome
  prot <- readRDS(file = paste0(path.data, "proteome_1", 
                                ifelse(time.point == 1, 
                                       "A", "B")))$data %>%
    dplyr::select(tidyselect::ends_with("_p")) %>% scale() %>%
    tibble::as_tibble()
  prot.feats <- readRDS(file = paste0(path.data, "proteome_1", 
                                      ifelse(time.point == 1, 
                                             "A", "B")))$feature.data
  internal.plotCorr(prot, prot.feats, 
                    paste0("corrPlot_prot", time.point))
  
  ## Methylome?
}

##### Function to load and filter fitted networks
clean.net <- function(path.corr, type.net, key.save = "") {
  path.save <- "results/images/"
  is.directed <- FALSE
  
  if (is.character(path.corr)) {
    conn <- gzfile(path.corr)
    df <- base::readRDS(conn) %>% tibble::as_tibble()
    net <- tidygraph::as_tbl_graph(df, directed = is.directed)
    close(conn)
  } else {
    net <- tidygraph::as_tbl_graph(path.corr, directed = is.directed)
  }
  
  net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(label = strsplit(name, "_")[[1]][2])
  
  net <- net %>%
    tidygraph::activate(., what = "nodes") %>%
    tidygraph::mutate(label = factor(net.nodes$label))
  
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
  ggplot2::ggsave(paste0(path.save, type.net, "_", key.save, ".png"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
}

##### Function to compute network properties
net.properties <- function(list.net.dfs, is.merged, how.to.join, 
                           path.save, key.save = "") {

  properties.all <- list()
  df.all <- list()
  for (n in names(list.net.dfs)) {
    l <- list.net.dfs[[n]]
    
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
      stop()
      df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                pcor, pval, qval, direction))
    } else {
      if (how.to.join == "inner") {
        colnames(df) <- c("from", "type.x", "idx.x", 
                          "to", "type.y", "idx.y", 
                          "type.x.2", "type.y.2", 
                          "pcor.x", "pval.x", "qval.x", 
                          "direction", 
                          "pcor.y", "pval.y", "qval.y")
        df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                  pcor.x, pcor.y, 
                                  pval.x, pval.y, 
                                  qval.x, qval.y, direction))
        return(df)
        stop()
      } else {
        colnames(df) <- c("from", "type.x", "idx.x", "to", "type.y", 
                          "idx.y", "type.x.2", "type.y.2", 
                          "pcor", "pval", "qval", "direction")
        stop()
        df <- dplyr::select(df, c(from, to, type.x, type.y, 
                                  pcor, pval, qval, direction))
      }
    }
    
    df$code <- as.character(n)
    if (dim.old != dim(df)[1]) {stop(call. = TRUE)}
    
    df.all <- append(df.all, list(df))
    
  } # End loop over networks
  
  # Create summary of properties
  df <- df.all %>% purrr::reduce(rbind)
  properties <- df %>%
    dplyr::select(-c(from, to)) %>%
    gtsummary::tbl_summary(by = "") %>%
    gtsummary::as_gt() %>%
    gt::gtsave(filename = paste0(path.save, "props_nets_", key.save, ".png"), 
               vwidth = 1500)
  
  return(properties)
} # End function compute properties networks

# properties.time.specific <- net.properties(list.net.dfs = list(
#   "visit A" = clean.net(path.corr = "../data/correlations/pcorr_mat_raw", 
#                         type.net = "time", key.save = "tmpA"), 
#   "visit B" = clean.net(path.corr = "../data/correlations/pcorr_mat_raw", 
#                         type.net = "time", key.save = "tmpB")
# ), is.merged = FALSE, how.to.join = "inner", 
# path.save = "results/images/", key.save = "time")

# properties.time.specific <- net.properties(list.net.dfs = list(
#   "merged" = clean.net(path.corr = "../data/correlations/merged_net", 
#                        type.net = "time", key.save = "tmp")
# ), is.merged = TRUE, how.to.join = "inner", 
# path.save = "results/images/", key.save = "merged")
