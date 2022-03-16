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

##### Tables #####
source("code/paper/tables.R")
## Creates table with population description
tables.population.desc()

## Describes networks
tables.networks(path.to.net = "../data/intermediate_res_ggm/merged_net.RData", 
                time.point = NULL, path.save = "results/images/")
tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 1, path.save = "results/images/")
tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 2, path.save = "results/images/")

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
