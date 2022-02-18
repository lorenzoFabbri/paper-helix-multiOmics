# Script to create .csv file containing the type of molecule for each
# molecule name (e.g., exposure -> "e"), and create tidied name

library(tidyverse)

source("code/multivariate_analysis/dictionaries.R")

path.meta <- "data/"
path.data <- "../data/"
path.data.methylome <- "../data/methylome/"
name.methylome <- "methylome_panel_ComBatSlide_6cells_v4.csv"
params <- list(
  time.point = 1
)

exposome <- readRDS(file = paste0(path.data, "exposome_1", 
                                  ifelse(params$time.point == 1, 
                                         "A", "B"))) %>%
  .$data %>% tibble::as_tibble() %>%
  dplyr::mutate(HelixID = from.sample.to.helix(SampleID)) %>%
  dplyr::select(tidyselect::ends_with("_e")) %>%
  colnames() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(label = "e") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(node.clean = strsplit(value, "_e")[[1]][1]) %>%
  dplyr::mutate(node.clean = strsplit(node.clean, "log.")[[1]][2]) %>%
  `colnames<-`(c("node", "label", "node.clean"))

prot <- readRDS(file = paste0(path.data, "proteome_1", 
                              ifelse(params$time.point == 1, 
                                     "A", "B"))) %>%
  .$data %>%
  dplyr::select(tidyselect::ends_with("_p")) %>%
  colnames() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(label = "p") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(node.clean = strsplit(value, "_p")[[1]][1]) %>%
  `colnames<-`(c("node", "label", "node.clean"))
metabs <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                ifelse(params$time.point == 1, 
                                       "A", "B"))) %>%
  .$data %>%
  dplyr::select(tidyselect::ends_with("_ms")) %>%
  colnames() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(label = "ms") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(node.clean = strsplit(value, "_ms")[[1]][1]) %>%
  `colnames<-`(c("node", "label", "node.clean"))
metabu <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                ifelse(params$time.point == 1, 
                                       "A", "B"))) %>%
  .$data %>%
  dplyr::select(tidyselect::ends_with("_mu")) %>%
  colnames() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(label = "mu") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(node.clean = strsplit(value, "_mu")[[1]][1]) %>%
  `colnames<-`(c("node", "label", "node.clean"))

# In order to avoid errors later on, I am renaming Creatinine and Taurine 
# since they are present both in the Serum and Urinary Metabolomes
metabu <- metabu %>%
  dplyr::mutate(node.clean = ifelse(node.clean == "Creatinine", "CreatinineU", node.clean)) %>%
  dplyr::mutate(node.clean = ifelse(node.clean == "Taurine", "TaurineU", node.clean))
metabs <- metabs %>%
  dplyr::mutate(node.clean = ifelse(node.clean == "Creatinine", "CreatinineS", node.clean)) %>%
  dplyr::mutate(node.clean = ifelse(node.clean == "Taurine", "TaurineS", node.clean))

if (length(intersect(metabs$node.clean, metabu$node.clean)) > 0) {stop(call. = TRUE)}

meth <- data.table::fread(paste0(path.data.methylome, name.methylome), 
                          nThread = 10) %>%
  dplyr::as_tibble() %>%
  dplyr::slice(1:3) %>%
  dplyr::select(tidyselect::starts_with("cg")) %>%
  dplyr::rename_with(., ~ paste0(., "_me") %>% unname() %>% unlist(), 
                     tidyselect::starts_with("cg")) %>%
  colnames() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(label = "me") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(node.clean = strsplit(value, "_me")[[1]][1]) %>%
  `colnames<-`(c("node", "label", "node.clean"))

label.mols <- dplyr::bind_rows(exposome, prot, metabs, metabu, meth)
label.mols %>%
  data.table::fwrite(., file = "data/labels_mols.csv", 
                     row.names = FALSE, col.names = TRUE)
if (dim(label.mols)[1] != (dim(exposome)[1]+dim(prot)[1]+dim(metabs)[1]+dim(metabu)[1]+dim(meth)[1])) {
  stop(call. = TRUE)
}
