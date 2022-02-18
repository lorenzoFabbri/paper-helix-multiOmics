# Script for generating images for paper

library(tidyverse)
library(tidygraph)

source("code/multivariate_analysis/dictionaries.R")

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
