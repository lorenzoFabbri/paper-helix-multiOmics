# Messy code to generate high-quality figures of the merged network
# and the various connected components

library(tidyverse)
library(igraph)
library(RColorBrewer)

rm(list=ls())
source("code/multivariate_analysis/dictionaries.R")

################################################################################
load("../data/intermediate_res_ggm/merged_net.RData")
net <- res$mod_2.2.2.5.5$net
num.ccs <- range(net %>% dplyr::as_tibble() %>% .$group)[2]

##### Add information on class of chemicals and -omic features
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
all.exposures <- dict.exposure.groups() %>%
  stack() %>% tibble::as_tibble() %>%
  `colnames<-`(c("exposure", "class")) %>%
  dplyr::mutate(class = as.character(class)) %>%
  dplyr::mutate(class = substr(class, 1, nchar(class) - 1))
all.exposures$exposure <- toupper(all.exposures$exposure)
net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
  tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.))) %>%
  dplyr::mutate(label = as.character(label)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(class = ifelse(
    label == "exposure", all.exposures[all.exposures$exposure == name, ]$class, 
    ifelse(
      label == "serum metabolome", metabs[metabs$Rvar == name, ]$Class, 
      ifelse(
        label == "proteome", "none", 
        "none"
      )
    )
  ))
#####

##### Database EWAS results
atlas <- data.table::fread("../data/methylome/results.txt")

path.res <- "results/focus_merged_ccs"
f <- list.files(path.res, include.dirs = T, full.names = T, recursive = T)
file.remove(f)
file.remove(f)

map.label.to.shape <- map.char.to.aes()[[1]]
map.class.to.col <- map.char.to.aes()[[2]]

##### Iterate over connected components of merged network
for (i in 1:num.ccs) {
  to.filter <- i
  
  load("../data/intermediate_res_ggm/merged_net.RData")
  net <- res$mod_2.2.2.5.5$net %>%
    tidygraph::mutate(class = as.factor(net.nodes$class)) %>%
    tidygraph::to_undirected()
  
  ##### Save to file to use with Cytoscape #####
  # Save seprarately edge list and node attributes
  if (!is.null(to.filter) & to.filter == 1) {
    # Save full network only once
    ig <- tidygraph::as.igraph(net)
    edge.list <- as.data.frame(igraph::get.edgelist(ig))
    vertex.attribs <- igraph::get.vertex.attribute(ig)
    vertex.attribs <- dplyr::bind_rows(
      c(as.data.frame(vertex.attribs$name), 
      as.data.frame(vertex.attribs$label), 
      as.data.frame(vertex.attribs$class))
    ) %>%
      `colnames<-`(c("node", "layer", "class"))
    write.csv(edge.list, 
              paste0(path.res, "/full_mergedNet", ".csv"), 
              row.names = FALSE, quote = FALSE)
    write.csv(vertex.attribs, 
              paste0(path.res, "/full_mergedNet_attributes", ".csv"), 
              row.names = FALSE, quote = FALSE)
  }
  
  if (!is.null(to.filter)) {
    # Filter specific connected component within the merged network
    net <- net %>% tidygraph::filter(group == to.filter)
  }
  
  # Not interested in pairwise connections alone
  if (dim(tidygraph::as_tibble(net))[1] < 3) {
    next
  }
  # Not interested in components with only Methylation sites
  if (length(unique(tidygraph::as_tibble(net)$label)) < 2) {
    if (as.character(unique(tidygraph::as_tibble(net)$label)) == "methylome") {
      next
    }
  }
  
  path.res.tmp <- paste0(path.res, "/cc_", to.filter)
  dir.create(path.res.tmp)
  ##### Save to file to use with Cytoscape #####
  # Save each connected component
  ig <- tidygraph::as.igraph(net)
  as.data.frame(igraph::get.edgelist(ig)) %>%
    write.csv(., paste0(path.res.tmp, "/cc_", to.filter, ".csv"), 
              row.names = FALSE, quote = FALSE)
  
  net <- net %>% tidygraph::activate("nodes") %>%
    # Group nodes based on the leading eigenvector of the modularity matrix
    tidygraph::mutate(gle = as.factor(tidygraph::group_leading_eigen(options = list(maxiter = 10000)))) %>%
    tidygraph::mutate(ce = as.numeric(tidygraph::centrality_eigen())) %>%
    tidygraph::mutate(degree = as.numeric(tidygraph::centrality_degree()))
  
  # Focus on each of the identified modules/clusters within the current cc
  num.gles <- net %>% tidygraph::as_tibble() %>% .$gle %>% as.numeric() %>% max()
  for (module in 1:num.gles) {
    net.module <- net %>% tidygraph::filter(gle == module) %>%
      tidygraph::mutate(gle.inner = as.factor(tidygraph::group_leading_eigen(options = list(maxiter = 10000))))
    
    # Tidy plot of module within current connected component
    gg.inner <- ggraph::ggraph(net.module, layout = "fr") +
      ggforce::geom_mark_hull(mapping = aes(x = x, y = y, 
                                            fill = class, 
                                            filter = !(class %in% c("none", 
                                                                    "phenols", 
                                                                    "pesticides", 
                                                                    "phthalates.high", 
                                                                    "phthalates.low"))), 
                              concavity = 1) +
      ggraph::geom_edge_link(colour = "darkgrey", alpha = 0.75) +
      ggraph::geom_node_point(mapping = aes(shape = label, size = degree, 
                                            color = class)) +
      ggraph::geom_node_text(mapping = aes(label = name), #filter = degree >= 5), 
                             repel = TRUE, size = 5) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_manual(values = map.class.to.col, drop = TRUE) +
      ggplot2::scale_shape_manual(values = map.label.to.shape, drop = TRUE) +
      ggplot2::theme(axis.line = ggplot2::element_blank(), 
                     axis.title = ggplot2::element_blank(), 
                     axis.text = ggplot2::element_blank(), 
                     axis.ticks = ggplot2::element_blank(), 
                     legend.title = element_text(size = 25), 
                     legend.text = element_text(size = 22)) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)), 
                      shape = ggplot2::guide_legend(override.aes = list(size = 3)), 
                      size = FALSE)
    ggplot2::ggsave(paste0(path.res.tmp, "/focus_cc_", to.filter, 
                           "mod_", module , ".pdf"), 
                    dpi = 720/2, 
                    width = 20, height = 12)
  } # End loop over modules within current connected component
  
  # Tidy plot of current connected component
  gg <- ggraph::ggraph(net, layout = "fr") +
    ggforce::geom_mark_hull(mapping = aes(x = x, y = y, 
                                          fill = class, 
                                          filter = !(class %in% c("none", 
                                                                  "phenols", 
                                                                  "pesticides", 
                                                                  "phthalates.high", 
                                                                  "phthalates.low"))), 
                            concavity = 1) +
    ggraph::geom_edge_link(colour = "darkgrey", alpha = 0.75) +
    ggraph::geom_node_point(mapping = aes(shape = label, size = degree, 
                                          color = class)) +
    ggraph::geom_node_text(mapping = aes(label = name, filter = degree >= 5), 
                           repel = TRUE, size = 5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = map.class.to.col, drop = TRUE) +
    ggplot2::scale_shape_manual(values = map.label.to.shape, drop = TRUE) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   axis.title = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), 
                   axis.ticks = ggplot2::element_blank(), 
                   legend.title = element_text(size = 25), 
                   legend.text = element_text(size = 22)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)), 
                    shape = ggplot2::guide_legend(override.aes = list(size = 3)), 
                    size = "none")
  ggplot2::ggsave(paste0(path.res.tmp, "/focus_cc_", to.filter, ".pdf"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
  
  # Find gene(s) associated to each CpG site
  cpgs <- net %>% dplyr::filter(label == "methylome") %>%
    as_tibble() %>% dplyr::select(name) %>% .$name %>% c()
  if (length(cpgs) > 0) {
    genes <- atlas %>% dplyr::filter(CpG %in% cpgs) %>%
      dplyr::select(Gene) %>% .$Gene %>%
      dplyr::as_tibble() %>%
      dplyr::filter(value != "-")
    if (dim(genes)[1] > 0) {
      genes %>% dplyr::mutate(value2 = strsplit(as.character(value), ";")) %>%
        tidyr::unnest(value2) %>% dplyr::select(value2) %>% .$value2 %>%
        dplyr::as_tibble() %>% dplyr::distinct() %>%
        write.csv(., paste0(path.res.tmp, "/genes_", to.filter, ".csv"), 
                  row.names = FALSE, quote = FALSE)
    }
  }
} # End loop over connected components of merged network
