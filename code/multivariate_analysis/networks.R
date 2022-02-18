# Author: Lorenzo Fabbri
# Script containing functions to generate and analyze networks

library(tidyverse)
library(mixOmics)
library(tidygraph)
library(ggraph)
library(gridExtra)
library(tidyHeatmap)
library(circlize)
library(viridis)

source("./code/ggm/ggm.R")

# Helper function to extract model.code (only the numbers)
extract.model.code <- function(string) {
  code <- strsplit(string, "_") %>% unlist(.) %>%
    .[2] %>% strsplit(., "") %>% unlist()
  
  return(code)
}

# Helper function to check threshold for constructing the networks
check.threshold.network <- function(model, comps) {
  cord.X <- cor(model$X, model$variates$X[, comps], use = "pairwise")
  cord.Y <- cor(model$Y, model$variates$X[, comps], use = "pairwise")
  
  ret <- cord.X %*% t(cord.Y) # Matrix of correlation values
  ret <- max(ret)
  return(ret)
}

# Function to generate/build networks
generate.networks <- function(model1, model2, params, type.network) {
  
  networks <- list()
  models <- list(model1, model2)
  
  if (type.network == "relevance") {
    dir.create(paste0("./results/networks/", type.network), 
               showWarnings = FALSE)
    
    # Relevance Networks from mixOmics
    for (mod in models) {
      comps <- 1:3
      threshold <- check.threshold.network(mod, comps)
      net.tmp <- mixOmics::network(mat = mod, comp = comps, 
                                   #cutoff = threshold - threshold / 4, 
                                   cutoff = 0, 
                                   show.edge.labels = TRUE, show.color.key = TRUE, 
                                   save = "jpeg", 
                                   name.save = paste0("./results/networks/", 
                                                      type.network, "/", "tmp"))
      networks <- append(networks, list(net.tmp))
    } # End loop models
  } # End Relevance Network
  
  return(networks)
} # End function generate.networks

# Helper function to create `tbl` graph object for later plotting
create.tbl <- function(net) {
  # Tidy labels exposures
  tmp <- tidygraph::as_tbl_graph(net$gR) %>% tibble::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(., labs = ifelse(stringr::str_detect(name, "log.hcp"), 
                                   unlist(strsplit(label, "log.hcp_,|_"))[2], 
                                   label))
  
  gg <- tidygraph::as_tbl_graph(net$gR) %>%
    tidygraph::activate(., what = "nodes") %>%
    dplyr::mutate(shape = ifelse(stringr::str_detect(name, "log.hcp"), 
                                 "exp", "omic")) %>%
    dplyr::mutate(labs = tmp$labs) %>%
    dplyr::mutate(community = as.factor(tidygraph::group_infomap())) %>%
    dplyr::mutate(centrality = igraph::centr_degree(.)$res) %>%
    tidygraph::activate(., what = "edges") %>%
    dplyr::mutate(sign = as.factor(sign(label)))
  
  return(gg)
}

# Helper function to create tidy plot of graph from `tidygraph` `tbl` object
plot.graph.tbl <- function(gg) {
  plt <- ggraph::ggraph(gg, layout = "stress") +
    ggraph::geom_edge_link(mapping = aes(width = pcor, 
                                         color = direction), 
                           alpha = 0.25) +
    ggraph::geom_node_point(mapping = aes(size = 3)) +
    ggraph::geom_node_text(mapping = aes(label = name), 
                           repel = TRUE) +
    ggraph::theme_graph()
  
  return(plt)
}

# Helper function to compare two networks (e.g., between time points)
compare.networks <- function(list.networks) {
  
  list.networks <- lapply(list.networks, create.tbl)
  net1 <- list.networks[[1]]
  net2 <- list.networks[[2]]
  
  # Create actual tibbles for ease of manipulation
  dfs <- list()
  for (net in list.networks) {
    df.tmp <- tidygraph::activate(net, what = "nodes") %>%
      igraph::as_data_frame() %>% tibble::as_tibble()
    communities <- tidygraph::activate(net, what = "nodes") %>%
      dplyr::select(community) %>% tibble::as_tibble()
    
    dfs <- append(dfs, list(df.tmp))
  }
  
  return(dfs)
} # End helper function compare networks

# Function to process generated networks
process.networks <- function(list.networks, params, type.network) {
  
  # Loop over list of networks (e.g., same networks but different time points)
  ret <- list()
  list.graphs <- list()
  for (net in list.networks) {
    gg <- create.tbl(net)
    list.graphs <- append(list.graphs, list(gg))
    
    # Plot graph
    plt <- plot.graph.tbl(gg)
    
    ret <- append(ret, list(plt))
  } # End loop over networks
  
  grid.plts <- gridExtra::arrangeGrob(grobs = ret, ncol = 2)
  ggplot2::ggsave(plot = grid.plts, 
                  filename = paste0(params$model.code, "_net", ".png"), 
                  device = "png", 
                  path = paste0("./results/networks/", type.network), 
                  dpi = 270, height = 12, width = 12)
  
} # End function process.networks

##################################################
##### FUNCTION TO PROCESS NETWORKS AS IN GGM #####
##################################################
process.networks.ggm <- function(list.networks, params) {
  
  nets.tidy <- list()
  i <- 1
  for (net in list.networks) {
    # Retrieve the matrix of "correlation" values
    pcor.mat <- net$M
    
    ## Tidy variables name
    # Exposures
    exp.groups <- dict.exposure.groups() %>% searchable::invert() %>% as.list()
    rownames(pcor.mat) <- rownames(pcor.mat) %>%
      lapply(., function(x) {
        if (startsWith(x, "log.hcp_")) {
          x <- strsplit(x, split = "_")[[1]][2]
          exp.group <- exp.groups[[x]] %>% as.character() %>%
            substr(., 1, nchar(.) - 1)
          paste0(x, "_", exp.group)
        } else { x }
      })
    
    # -omics
    omic.type <- params[[i]]$omic.type %>% as.character()
    colnames(pcor.mat) <- colnames(pcor.mat) %>%
      lapply(., function(x) {
        if (grepl("log.", x)) {
          strsplit(x, split = "log.")[[1]][2] %>%
            paste0(., "_", omic.type)
        } else {
          paste0(x, "_", omic.type)
        }
      })
    
    nets.tidy[[params[[i]]$model.code]] <- list(
      model = pcor.mat, 
      params = "", 
      sample.size = ""
    )
    
    i <- i + 1
  } # End loop over networks
  
  ## Process tidied networks
  params <- list(
    method.ggm = "qval", 
    cutoff.ggm = 0.05, 
    method.mixed = tidygraph::group_infomap, 
    path.save.plots = "./results/networks/relevance/"
  )
  
  # We use the same pipeline used for the GGMs, but these are 
  # *relevance networks*
  # The list of networks contains two networks: one for each time point
  processed.nets <- process.ggms(list.ggms = nets.tidy, active = "relevance", 
                                 filter.mixed.interactions = FALSE, 
                                 is.directed = FALSE, 
                                 params = params)
  
} # End function process.networks.ggm

################################################
##### Function to compute properties of RN #####
################################################
process.rn <- function(nets, to_keep) {
  
  # Iterate over fitted models
  nets.all <- list()
  for (n in names(nets)) {
    l <- nets[[n]]$net
    
    # code: scaling.type adjustment.type model.type time.point omic.type exposure.type
    code <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
      as.list() %>% .[[1]] %>% strsplit(., "") %>% .[[1]] %>% .[5:6] %>%
      paste0(., collapse = "")
    
    if (code != to_keep) next
    cat("\n", paste0("Processing GGM: ", n, "..."), "\n")
    
    nets.all[[n]]$net <- l
  } # End loop over GGMs
  
  net.props <- net.properties(ggms = nets.all, 
                              type.networks = "relevance", 
                              path.save = "results/networks/relevance/", 
                              code.save = to_keep, is.merged = FALSE)
  
}

# Function to run analyses
main <- function(path.models, type.network) {
  
  # Iterate over all directories containing saved images
  file.names <- list.files(path.models, 
                           full.names = TRUE)
  models1 <- list()
  models2 <- list()
  idx.time.point <- 4
  for (dir in file.names) {
    code <- strsplit(dir, "/")[[1]][6] %>% extract.model.code(.)
    
    if (code[idx.time.point] == "1") {
      models1 <- append(models1, dir)
    } else if (code[idx.time.point] == "2") {
      models2 <- append(models2, dir)
    } else {
      stop(call. = TRUE)
    }
  } # End loop over list directories
  
  processed.ggms.all <- list()
  merged.all <- list()
  for (idx in 1:length(models1)) {
    model1 <- get(load(paste0(models1[idx], "/ret.Rdata")))$model
    model2 <- get(load(paste0(models2[idx], "/ret.Rdata")))$model
    params1 <- get(load(paste0(models1[idx], "/pls.params.Rdata")))
    params2 <- get(load(paste0(models2[idx], "/pls.params.Rdata")))
    params <- list(params1, params2)
    
    # Run function to create networks
    networks <- generate.networks(model1, model2, params, type.network)
    
    # Run function to process networks
    #process.networks(networks, params, type.network)
    processed.ggms <- process.networks.ggm(list.networks = networks, 
                                           params = params)
    processed.ggms.all <- append(processed.ggms.all, processed.ggms)
    
    # Run function to merge networks by time point
    n <- names(processed.ggms)[1]
    code.omic <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
      as.list() %>% .[[1]] %>% strsplit(., "") %>% .[[1]] %>% .[5]
    code.exp <- strsplit(n, "_")[[1]][2] %>% strsplit(., "\\.") %>% .[[1]] %>%
      as.list() %>% .[[1]] %>% strsplit(., "") %>% .[[1]] %>% .[6]
    #code <- paste0("o", code.omic, "e", code.exp)
    res <- merge.networks(ggms = processed.ggms, 
                          exposure.group = code.exp, omic.type = code.omic, 
                          how.to.join = "inner", type.networks = "relevance", 
                          path.save = "./results/networks/relevance/")
    merged.all <- append(merged.all, res)
  } # End loop over models
  
  return(
    list(
      processed.ggms = processed.ggms.all, 
      merged.nets = merged.all
    )
  )
}
