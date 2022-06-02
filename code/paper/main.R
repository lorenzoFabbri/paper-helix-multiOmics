# Main script to run analyses and analyze results from cluster

rm(list=ls())

##### Figures #####
source("code/paper/images.R")
## QQ-plot of p-values from results EWAS
qq.pval.ewas(path = "results/ewaff_win_sva/")

## Visualize significant CpG sites by chemical class
visualize.common.cpgs(path = "results/ewaff_win_sva/", 
                      threshold = 0.2)

## Visualize results of bootstrapping merged GGM
plot.bootstrapping.nets(path = "results/ggm/paper_boot/", 
                        path.save = "results/images/")

## Visualize correlation exposome between time points
plot.cor.exp.time()

## Heatmaps of exposome and -omics
plot.heatmaps(time.point = 1)
plot.heatmaps(time.point = 2)

## Tidy and filter networks
### Full, merged network
clean.net(path.corr = "../data/correlations/merged_net", 
          type.net = "merged", key.save = "paper")

### Networks centered on each exposure
all.exposures <- dict.exposure.groups() %>%
  stack() %>% tibble::as_tibble() %>%
  `colnames<-`(c("exposure", "class")) %>%
  dplyr::mutate(class = as.character(class)) %>%
  dplyr::mutate(class = substr(class, 1, nchar(class) - 1)) %>%
  dplyr::mutate(exposure = toupper(exposure))
for (chem in all.exposures$exposure) {
  clean.net(path.corr = "../data/correlations/merged_net", 
            type.net = "merged", key.save = "paper", 
            center.node = chem)
}

## Visualize largest "component" of merged network with heatmap
load("../data/intermediate_res_ggm/merged_net.RData")
plot.net.as.heatmap(res$mod_2.2.2.5.5$net)

## Visualize results pathway enrichment from MetaboAnalyst
df <- read.csv("results/pathways/metaboanalyst_res/pathway_results.csv") %>%
  tibble::as_tibble()
plt <- ggplot2::ggplot(data = df, mapping = aes(x = Impact, 
                                                y = X.log10.p., 
                                                label = X)) +
  ggplot2::geom_point(mapping = aes(size = Impact, 
                                    colour = -X.log10.p.), 
                      show.legend = FALSE) +
  ggplot2::geom_text(check_overlap = TRUE, 
                     data = subset(df, X.log10.p. >= 5 | 
                                     Impact >= 0.5), 
                     hjust = "inward") +
  ggplot2::labs(x = "Pathway impact", 
                y = "-log10(p)") +
  ggplot2::scale_color_distiller(palette = "YlOrRd") +
  ggplot2::theme_minimal()
ggplot2::ggsave(filename = "results/images/pathway_enrichment.png", 
                height = 15, width = 15, dpi = 720, plot = plt)

## Distribution of partial correlations and histogram of p-values for
## processed networks
unproc.net1 <- readRDS("../data/correlations/network_raw_mod_2.2.1.5.5")
unproc.net2 <- readRDS("../data/correlations/network_raw_mod_2.2.2.5.5")
unproc.nets <- list(unproc.net1, unproc.net2) %>%
  dplyr::bind_rows(.id = "visit") %>%
  dplyr::mutate(visit = ifelse(visit == 1, "A", "B"))

tmp <- unproc.nets %>%
  ggplot2::ggplot(mapping = aes(x = pcor, 
                                fill = visit)) +
  ggplot2::geom_histogram(alpha = 0.4, position = "identity", 
                          show.legend = TRUE) +
  ggplot2::facet_wrap(~ visit) +
  ggplot2::theme_light() +
  ggplot2::theme(strip.background = element_blank(), 
                 panel.border = element_blank(), 
                 panel.grid = element_blank(), 
                 axis.line = element_line(), 
                 text = ggplot2::element_text(size = 20)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(x = "partial correlations") +
  ggplot2::scale_fill_brewer(palette = "Accent")
tmp.pvals <- unproc.nets %>%
  ggplot2::ggplot(mapping = aes(x = pval, 
                                fill = visit)) +
  ggplot2::geom_histogram(alpha = 0.4, position = "identity", 
                          show.legend = FALSE) +
  ggplot2::facet_wrap(~ visit) +
  ggplot2::theme_light() +
  ggplot2::theme(strip.background = element_blank(), 
                 panel.border = element_blank(), 
                 panel.grid = element_blank(), 
                 axis.line = element_line(), 
                 text = ggplot2::element_text(size = 20)) +
  ggplot2::labs(x = "p-values") +
  ggplot2::scale_fill_brewer(palette = "Accent")
plt <- egg::ggarrange(tmp.pvals, tmp, ncol = 1, heights = c(0.3, 1.2))
ggplot2::ggsave(filename = "results/final_material_paper_v2/pcorsPvals.png", 
                height = 10, width = 15, dpi = 720, plot = plt)

## Distribution of types of edges by chemical class for merged network
all.exposures <- dict.exposure.groups() %>%
  stack() %>% tibble::as_tibble() %>%
  `colnames<-`(c("exposure", "class")) %>%
  dplyr::mutate(class = as.character(class)) %>%
  dplyr::mutate(class = substr(class, 1, nchar(class) - 1)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(class = ifelse(
    class %in% c("phthalates.low", "phthalates.high"), 
    "phthalates", class
  )) %>% dplyr::mutate(exposure = toupper(exposure))
merged.net <- get(load("../data/intermediate_res_ggm/merged_net.RData"))
merged.net <- tidy.graph(merged.net$mod_2.2.2.5.5$net)
merged.net <- dplyr::left_join(merged.net, all.exposures, 
                               by = c("name.x" = "exposure"))
merged.net <- dplyr::left_join(merged.net, all.exposures, 
                               by = c("name.y" = "exposure"))
merged.net <- merged.net %>%
  dplyr::mutate(chem.class = dplyr::coalesce(class.x, class.y)) %>%
  dplyr::select(-c(class.x, class.y)) %>%
  dplyr::mutate(chem.class = as.factor(chem.class))
list.pies <- list()
for (chem in levels(merged.net$chem.class)) {
  df.tmp <- merged.net %>%
    dplyr::filter(chem.class %in% c(NA, chem)) %>%
    dplyr::select(name.x, name.y, label.x, label.y)
  df.tmp <- tibble::tibble(node = c(df.tmp$name.x, df.tmp$name.y), 
                           label = c(df.tmp$label.x, df.tmp$label.y)) %>%
    dplyr::filter(label != "exposure") %>%
    dplyr::distinct() %>%
    dplyr::group_by(label) %>% dplyr::summarise(freq = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(perc = round(freq / sum(freq), 2)) %>%
    dplyr::mutate(perc = scales::percent(perc))
  pie <- ggplot2::ggplot(df.tmp, mapping = aes(x = "", 
                                               y = freq, fill = label)) +
    ggplot2::geom_bar(width = 1, stat = "identity") +
    ggplot2::geom_text(aes(label = perc), 
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::theme_void() + ggplot2::scale_fill_brewer()
  list.pies <- append(list.pies, list(pie))
}

## Code to convert Table 3 to figure(s)
counts <- read.csv("results/final_material_paper_v2/edge_countTable.csv") %>%
  tibble::as_tibble() %>%
  dplyr::select(layer.node.x, 
                layer.node.y, 
                n, n.A, n.B, avg.pc) %>%
  tidyr::pivot_longer(cols = starts_with("n"), 
                      names_to = "n", values_to = "count") %>%
  dplyr::mutate(layer.node.x = dplyr::case_when(
    layer.node.x == "exposure" ~ "e", 
    layer.node.x == "serum metabolite" ~ "sm", 
    layer.node.x == "urinary metabolite" ~ "um", 
    layer.node.x == "protein" ~ "p", 
    layer.node.x == "methylome" ~ "m", 
  )) %>%
  dplyr::mutate(layer.node.y = dplyr::case_when(
    layer.node.y == "exposure" ~ "e", 
    layer.node.y == "serum metabolite" ~ "sm", 
    layer.node.y == "urinary metabolite" ~ "um", 
    layer.node.y == "protein" ~ "p", 
    layer.node.y == "methylome" ~ "m", 
  )) %>%
  tidyr::unite("edge_type", layer.node.x:layer.node.y, 
               remove = TRUE, sep = "-") %>%
  dplyr::mutate(n = ifelse(n == "n", "merged", 
                           ifelse(n == "n.A", "visit A", "visit B"))) %>%
  dplyr::mutate(avg.pc = paste0(avg.pc, "%"))
plt <- ggplot2::ggplot(data = counts, mapping = aes(x = edge_type, 
                                                    y = count, 
                                                    fill = n)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::theme_light() +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.5), 
                 text = ggplot2::element_text(size = 20)) +
  ggplot2::labs(x = "", y = "counts", fill = "network") +
  ggplot2::geom_text(aes(label = ifelse(n == "merged", avg.pc, "")), 
                     position = ggplot2::position_dodge(width = 0.9), 
                     color = "black", vjust = -0.3) +
  ggplot2::scale_fill_brewer(palette = "Accent")
ggplot2::ggsave(filename = "results/final_material_paper_v2/barplots_countsEdges.png", 
                height = 10, width = 15, dpi = 720, plot = plt)

################################################################################
##### Tables #####
################################################################################
source("code/paper/tables.R")
## Creates table with population description
tables.population.desc()

## Describes networks
merged <- tables.networks(path.to.net = "../data/intermediate_res_ggm/merged_net.RData", 
                time.point = NULL, path.save = "results/images/")
net1 <- tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 1, path.save = "results/images/")
net2 <- tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 2, path.save = "results/images/")
# Table describing networks for T1, T2 and merged (no sub-graphs for chemicals)
gtsummary::tbl_merge(list(net1, net2), 
                     tab_spanner = c("A", "B")) %>%
  gtsummary::as_gt() %>%
  gt::gtsave(file = paste0("results/images/desc_allTimesPaper.rtf"))
merged %>%
  gtsummary::as_gt() %>%
  gt::gtsave(file = paste0("results/images/desc_mergedPaper.rtf"))

# Compare top, mixed associations across time points (from time-specific 
# networks)
assoc1 <- read.csv("results/images/mixedEdges_net1.csv") %>% tibble::as_tibble() %>%
  dplyr::mutate(visit = "A")
assoc2 <- read.csv("results/images/mixedEdges_net2.csv") %>% tibble::as_tibble() %>%
  dplyr::mutate(visit = "B")
N = 5
dplyr::bind_rows(assoc1[1:N, ], assoc2[1:N, ]) %>%
  write.csv("results/images_locked/top5Mixed_timeNets.csv", 
            quote = FALSE, row.names = FALSE)
assoc.all.a <- dplyr::full_join(assoc1, assoc2, by = c("node.a" = "node.a", 
                                                      "node.b" = "node.b", 
                                                      "layer.a" = "layer.a", 
                                                      "layer.b" = "layer.b"))
assoc.all.b <- dplyr::full_join(assoc1, assoc2, by = c("node.a" = "node.b", 
                                                      "node.b" = "node.a", 
                                                      "layer.a" = "layer.b", 
                                                      "layer.b" = "layer.a"))
assoc.all <- dplyr::bind_rows(assoc.all.a, assoc.all.b)

## Tables summarizing difference among time-specific and merged networks
merged.net <- get(load("../data/intermediate_res_ggm/merged_net.RData"))
merged.net <- merged.net$mod_2.2.2.5.5$net
processed.nets <- get(load("../data/intermediate_res_ggm/processed_ggms.RData"))
t1.net <- processed.nets$mod_2.2.1.5.5$net
t2.net <- processed.nets$mod_2.2.2.5.5$net
summary.net <- function(net) {
  num.nodes <- igraph::vcount(net)
  num.edges <- igraph::ecount(net)
  diam <- igraph::diameter(net, directed = FALSE)
  edge.density <- igraph::edge_density(net)
  
  return(list(
    num.nodes = num.nodes, 
    num.edges = num.edges, 
    diameter = diam, 
    edge.density = round(edge.density, 2)
  ))
}
s.me <- summary.net(merged.net)
s.t1 <- summary.net(t1.net)
s.t2 <- summary.net(t2.net)
summary.nets <- list(merged = tibble::as_tibble(s.me), 
                     T1 = tibble::as_tibble(s.t1), 
                     T2 = tibble::as_tibble(s.t2)) %>%
  purrr::reduce(dplyr::bind_rows)
colnames(summary.nets) <- c("number nodes", "number edges", 
                            "diameter", "edge density")
summary.nets[["network"]] <- c("merged", "T1", "T2")
summary.nets <- dplyr::relocate(summary.nets, "network") %>%
  gt::gt() %>%
  gt::gtsave(filename = "results/images/comparison_nets.rtf")

##### Convert csv files for paper to rtf files to be added to all_tables file
read.csv("results/final_material_paper/edge_countTable.csv") %>%
  tibble::as_tibble() %>%
  `colnames<-`(c("layer node x", "layer node y", "n", "p", 
                 "n A", "p A", "n B", "p B", "average change (%)")) %>%
  gt::gt() %>%
  gt::tab_style(gt::cell_text(weight = "bold"), 
                locations = gt::cells_column_labels()) %>%
  gt::cols_align(align = "center") %>%
  gt::gtsave(., filename = "results/final_material_paper/edge_countTable.rtf", 
             zoom = 4, expand = 7)

read.csv("results/final_material_paper/top5Mixed_timeNets.csv", 
         header = TRUE, check.names = FALSE) %>%
  tibble::as_tibble() %>%
  `colnames<-`(c("node x", "layer x", "degree x", 
                 "node y", "layer y", "degree y", 
                 "pcor", "qval", "visit")) %>%
  gt::gt() %>%
  gt::tab_style(gt::cell_text(weight = "bold"), 
                locations = gt::cells_column_labels()) %>%
  gt::cols_align(align = "center") %>%
  gt::gtsave(., filename = "results/final_material_paper/top5Mixed_timeNets.rtf", 
             zoom = 4, expand = 7)

read.csv("results/final_material_paper/mixedEdges_netMerged_locked.csv", 
         header = TRUE, check.names = FALSE) %>%
  tibble::as_tibble() %>%
  `colnames<-`(c("chemical class", "chemical", "omic layer", 
                 "omic feature", "pcor A", "pcor B", 
                 "qval A", "qval B", "genes")) %>%
  gt::gt() %>%
  gt::tab_style(gt::cell_text(weight = "bold"), 
                locations = gt::cells_column_labels()) %>%
  gt::cols_align(align = "center") %>%
  gt::gtsave(., filename = "results/final_material_paper/mixedEdges_netMerged_locked.rtf", 
             zoom = 4, expand = 7)
