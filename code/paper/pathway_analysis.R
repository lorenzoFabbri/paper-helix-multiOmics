# Functions to generate data as input for MetaboAnalyst to perform
# pathway analysis

rm(list = ls())
library(tidyverse)

source("code/multivariate_analysis/dictionaries.R")

####################################################################
##### Generate list of metabolites for specific chemical class #####
####################################################################
generate.list.mets <- function(chemical.class, time.point) {
  
  path.data <- "../data/"
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                  ifelse(time.point == 1, 
                                         "A", "B"))) %>% .$feature.data %>%
    dplyr::select(c(Rvar, var, CHEBI))
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                  ifelse(time.point == 1, 
                                         "A", "B"))) %>% .$feature.data %>%
    dplyr::select(c(Rvar, var, CHEBI))
  
  # Create vector of all exposures
  all.exposures <- dict.exposure.groups() %>%
    stack() %>% tibble::as_tibble() %>%
    `colnames<-`(c("exposure", "class")) %>%
    dplyr::mutate(class = as.character(class)) %>%
    dplyr::mutate(class = substr(class, 1, nchar(class) - 1)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(class = ifelse(
      class %in% c("phthalates.low", "phthalates.high"), 
      "phthalates", class
    ))
  
  # Load time-specific network to select metabolites associated to
  # chemicals of interest
  net <- get(load("../data/intermediate_res_ggm/processed_ggms.RData"))
  net <- net[[time.point]]$net
  nodes <- net %>% tidygraph::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(name = ifelse(
      type == "e", 
      strsplit(name, c("log.|_"))[[1]][2], 
      name
    ))
  nodes <- dplyr::full_join(nodes, all.exposures, by = c("name" = "exposure"))
  net <- net %>%
    tidygraph::mutate(class = nodes$class) %>%
    tidygraph::mutate(class = tidyr::replace_na(class, "omic"))
  
  # Filter exposures based on given chemical class
  exposures <- all.exposures %>%
    dplyr::filter(class == chemical.class) %>%
    dplyr::select(exposure) %>%
    c() %>% unname() %>% unlist()
  
  df <- tidy.graph(net %>% tidygraph::filter(class %in% c(chemical.class, 
                                                          "omic"))) %>%
    dplyr::select(c(name.x, name.y))
  
  list.metabolites <- c(df$name.x, df$name.y)
  if (length(list.metabolites) != nrow(df) * 2) {stop(call. = TRUE)}
  list.metabolites <- list.metabolites %>%
    tibble::as_tibble() %>% `colnames<-`(c("Rvar.dirty")) %>%
    dplyr::filter(grepl("_mu|_ms", Rvar.dirty)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Rvar = strsplit(Rvar.dirty, "_")[[1]][1]) %>%
    dplyr::distinct()
  list.metabolites <- dplyr::left_join(list.metabolites, metabs) %>%
    dplyr::mutate(CHEBI = as.character(CHEBI))
  list.metabolites <- dplyr::left_join(list.metabolites, metabu, by = "Rvar")
  
  list.metabolites.clean <- list.metabolites %>%
    dplyr::mutate(var = dplyr::coalesce(var.x, var.y), .keep = "unused") %>%
    dplyr::mutate(CHEBI = dplyr::coalesce(CHEBI.x, CHEBI.y), .keep = "unused")
  # Convert CHEBI to CID
  # dplyr::mutate(pcid = webchem::cts_convert(query = paste0("CHEBI:", CHEBI), 
  #                                           from = "chebi", to = "pubchem cid") %>%
  #                 unname() %>% .[[1]] %>% .[1]) %>%
  # Convert from CHEBI to KEGG
  # dplyr::mutate(kegg = webchem::cts_convert(query = paste0("CHEBI:", CHEBI), 
  #                                           from = "chebi", to = "kegg") %>%
  #                 unname() %>% .[[1]] %>% .[1])
  
  list.for.metaboanalyst <- list.metabolites.clean %>%
    #tidyr::unnest(chemical.name) %>%
    dplyr::distinct(var, .keep_all = TRUE) %>%
    tidyr::drop_na(var) %>% dplyr::select(var)
  
  # Merge with manually-created txt file of "expanded" chemical names (for those
  # not founded by MetaboAnalyst)
  mapping <- read.table("results/pathways/mapping_metabolites.txt", 
                        header = TRUE, sep = ";") %>%
    tibble::as_tibble()
  list.for.metaboanalyst <- list.for.metaboanalyst %>%
    dplyr::full_join(mapping, by = c("var" = "old"))
  list.for.metaboanalyst$new <- ifelse(
    is.na(list.for.metaboanalyst$new), 
    list.for.metaboanalyst$var, list.for.metaboanalyst$new
  )
  list.for.metaboanalyst %>%
    .$new %>% c() %>% noquote() %>% cat(., sep = "\n", 
                                        file = paste0("results/pathways/", 
                                                      "time", time.point, 
                                                      "class_", chemical.class, 
                                                      ".txt"))
}

chemicals <- c("phenols", "phthalates", "pesticides")
time.points <- c(1, 2)
#####
#generate.list.mets(chemical.class = chemicals[1], time.point = time.points[1])
#####
for (time.point in time.points) {
  for (chemical in chemicals) {
    cat("\n=======================================================\n")
    cat(paste0("Time point: ", time.point, " - Chemical class: ", 
               chemical, ".\n"))
    generate.list.mets(chemical.class = chemical, 
                       time.point = time.point)
    cat("=======================================================\n")
    #readline(prompt = "Press [enter] to continue...")
  }
}
