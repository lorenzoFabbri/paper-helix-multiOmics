# Main script to run analyses and analyze results from cluster

rm(list=ls())
source("code/paper/images.R")
source("code/paper/tables.R")

##### Figures #####
## QQ-plot of p-values from results EWAS
qq.pval.ewas(path = "results/ewaff_win_sva/")

## Visualize significant CpG sites by chemical class
visualize.common.cpgs(path = "results/ewaff_win_sva/", 
                      threshold = 0.2)

## Visualize results of bootstrapping merged GGM
plot.bootstrapping.nets(path = "results/ggm/paper_boot/")

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

##### Tables #####
## Creates table with population description
tables.population.desc()

## Describes networks
tables.networks(path.to.net = "../data/intermediate_res_ggm/merged_net.RData", 
                time.point = NULL)
tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 1)
tables.networks(path.to.net = "../data/intermediate_res_ggm/processed_ggms.RData", 
                time.point = 2)
