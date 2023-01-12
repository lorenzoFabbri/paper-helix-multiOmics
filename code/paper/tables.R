# Script for generating tables for paper

library(tidyverse)

source("code/multivariate_analysis/dictionaries.R")

##### Helper function to load metadata by time point
load.metadata <- function(time.point) {
  path.meta <- "data/"
  
  # For the first review of Env Int, I was asked to add info
  # on smoking and maternal education. I found these variables
  # in a dataset from another project in the HPC cluster
  .file.name <- "HELIX_subcohort_raw.csv"
  extra.covars <- readr::read_csv(paste0("../../../analyses/GE_BMI_IF/db/", 
                                  .file.name), 
                                  col_names = TRUE, 
                                  col_types = readr::cols()) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::select(-tidyselect::contains("HelixID")) %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID)) %>%
    dplyr::select(c(HelixID, 
                    hs_globalexp_None, hs_globalexp2_None, 
                    hs_smk_parents_None, hs_globalsmok_m_None, 
                    h_edumc_None))
  
  season <- readr::read_csv(paste0(path.meta, "metadata_old.csv"), 
                            col_names = TRUE, 
                            col_types = readr::cols()) %>%
    dplyr::select(c(SampleID, season, period)) %>%
    dplyr::filter(period == ifelse(time.point == 1, "A", "B")) %>%
    dplyr::mutate(SampleID = gsub("EDP", "EDE", SampleID)) %>%
    dplyr::select(-c(period)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(HelixID = substr(SampleID, 1, nchar(SampleID) - 3)) %>%
    dplyr::select(-c(SampleID))
  
  metadata <- readr::read_csv(file = paste0(path.meta, "meta", 
                                            time.point, ".csv"), 
                              col_names = TRUE, 
                              col_types = readr::cols()) %>%
    dplyr::select(-tidyselect::contains("HelixID")) %>%
    dplyr::mutate(HelixID = from.sample.to.helix(SampleID)) %>%
    dplyr::inner_join(season)
  old.nrows <- nrow(metadata)
  metadata <- metadata %>%
    dplyr::inner_join(extra.covars, by = "HelixID")
  if (old.nrows != nrow(metadata)) { stop("Number of rows is different.", 
                                          call. = TRUE) }
  
  common.subjects <- readr::read_csv(paste0(path.meta, 
                                            paste0("common_samples_t", 
                                                   time.point, ".csv")), 
                                     col_types = readr::cols()) %>%
    `colnames<-`(c("HelixID"))
  
  # Final dataset
  metadata <- metadata %>%
    dplyr::filter(HelixID %in% common.subjects$HelixID) %>%
    dplyr::rename(zBMI = hs_zbmi_theano) %>%
    dplyr::select(-c("SampleID"))
  
  return(metadata)
}

##### Function to create table with population description
tables.population.desc <- function(to.save, tbl.by = "Period") {
  meta.1 <- load.metadata(time.point = 1)
  meta.2 <- load.metadata(time.point = 2)
  metadata <- dplyr::bind_rows(meta.1, meta.2)
  
  # Add information regarding difference in months between visits
  merged <- dplyr::inner_join(meta.1, meta.2, by = "HelixID")
  merged <- merged %>%
    dplyr::mutate(diff.visits = age_sample_months.x.y - age_sample_months.x.x)
  metadata <- dplyr::full_join(metadata, 
                               merged %>%
                                 dplyr::select(HelixID, diff.visits), 
                               by = "HelixID")
  metadata <- metadata %>%
    dplyr::mutate(season = factor(season) %>% forcats::fct_explicit_na(), 
                  h_edumc_None = factor(h_edumc_None) %>%
                    forcats::fct_explicit_na(), 
                  hs_globalsmok_m_None = factor(hs_globalsmok_m_None) %>%
                    forcats::fct_explicit_na()) %>%
    dplyr::mutate(h_edumc_None = forcats::fct_recode(
      h_edumc_None, 
      "Low" = "1", "Middle" = "2", "High" = "3"
    ))
  
  # Table
  hs_mealblood <- "Time to last meal (hours)"
  diff.visits <- "Difference between visits (months)"
  h_edumc_None <- "Maternal education"
  hs_globalsmok_m_None <- "Tobacco exposure"
  # Helper function to generate summary table
  .summary <- function(metadata) {
    gtsummary::tbl_summary(metadata %>%
                             dplyr::select(-c(cohort)) %>%
                             dplyr::rename(., Cohort = cohort.x) %>%
                             dplyr::rename(., Sex = e3_sex.x) %>%
                             dplyr::rename(., Age = age_sample_years.x) %>%
                             dplyr::rename(., Season = season) %>%
                             dplyr::rename(., Ethnicity = h_ethnicity_c.x) %>%
                             dplyr::mutate(., Period = replace(
                               Period, Period == "1A", "A"
                             )) %>%
                             dplyr::mutate(., Period = replace(
                               Period, Period == "1B", "B"
                             )) %>%
                             dplyr::mutate(Ethnicity = factor(
                               Ethnicity, 
                               levels = c("Caucasian", "Pakistani", "Other")
                             ), 
                             Period = as.factor(Period)), 
                           include = c("Cohort", "Sex", "Age", "diff.visits", 
                                       "Ethnicity", "zBMI", 
                                       "hs_dift_mealblood_imp", "Season", 
                                       "h_edumc_None", "hs_globalsmok_m_None"
                           ), 
                           label = list(
                             hs_dift_mealblood_imp ~ hs_mealblood, 
                             diff.visits ~ diff.visits, 
                             h_edumc_None ~ h_edumc_None, 
                             hs_globalsmok_m_None ~ hs_globalsmok_m_None
                           ), missing = "ifany", 
                           by = tbl.by) %>%
      gtsummary::bold_labels()
  }
  if (tbl.by == "Cohort") {
    .path <- "../../lorenzoF_phd/PhD/papers/paper1_helix_multiOmics/EnvInt_Sep2022/review1/dat/"
    infos <- lapply(metadata %>% dplyr::group_split(Period), 
                    FUN = .summary)
    if (to.save == TRUE) {
      gt::gtsave(infos[[1]] %>% gtsummary::as_gt(), 
                 filename = paste0(.path, "pop_desc_a.rtf"), 
                 zoom = 4, expand = 7)
      gt::gtsave(infos[[2]] %>% gtsummary::as_gt(), 
                 filename = paste0(.path, "pop_desc_b.rtf"), 
                 zoom = 4, expand = 7)
    } else { return(infos) }
  } else {
    info <- .summary(metadata = metadata)
    if (to.save == TRUE) {
      gt::gtsave(info %>% gtsummary::as_gt(), 
                 filename = paste0(.path, "pop_desc.rtf"), 
                 zoom = 4, expand = 7)
    } else { return(info) }
  }
}

##### Helper function for `tables.networks` to process networks
.tables.networks <- function(net, time.point, 
                             chem.classes, chem.class) {
  if (!is.null(chem.class)) {
    net.sub <- net %>%
      tidygraph::filter(!(class %in% setdiff(chem.classes, chem.class)))
  } else {
    net.sub <- net
  }
  net.sub.old <- net.sub
  net.sub <- tidy.graph(net.sub)
  
  df.filtered <- net.sub %>%
    dplyr::select(-dplyr::any_of(c("class.x", "class.y", 
                                   "pval.x", "pval.y", "pval", "prob"))) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), round, 3))
  
  if (is.null(time.point)) {
    # Merged network
    df.filtered <- df.filtered %>%
      dplyr::arrange(name.x, name.y, label.x, label.y, 
                     desc(pcor.x), desc(pcor.y)) %>%
      `colnames<-`(c("node.a", "layer.a", "degree.a", 
                     "node.b", "layer.b", "degree.b", 
                     "pcor.a", "qval.a", "sign", "pcor.b", "qval.b")) %>%
      dplyr::select(c(node.a, layer.a, degree.a, 
                      node.b, layer.b, degree.b, 
                      sign, 
                      pcor.a, qval.a, 
                      pcor.b, qval.b))
    colnames(df.filtered) <- c("node.a", "layer.a", "degree.a", 
                               "node.b", "layer.b", "degree.b", 
                               "sign", 
                               "pcor.t1", "qval.t1", "pcor.t2", "qval.t2")
  } else {
    # Time-specific network
    df.filtered <- df.filtered %>%
      dplyr::arrange(name.x, name.y, label.x, label.y, 
                     desc(pcor)) %>%
      `colnames<-`(c("node.a", "layer.a", "degree.a", 
                     "node.b", "layer.b", "degree.b", 
                     "pcor", "qval", "sign")) %>%
      dplyr::select(c(node.a, layer.a, degree.a, 
                      node.b, layer.b, degree.b, 
                      sign, 
                      pcor, qval))
  }
  
  # Proportion of layer types by node
  n.a <- df.filtered %>% dplyr::select(c(node.a, layer.a)) %>%
    `colnames<-`(c("node", "layer"))
  n.b <- df.filtered %>% dplyr::select(c(node.b, layer.b)) %>%
    `colnames<-`(c("node", "layer"))
  n.all <- dplyr::bind_rows(n.a, n.b) %>%
    dplyr::distinct()
  
  # Add infomation on degree of nodes
  old.dim <- nrow(n.all)
  degree.nodes <- net.sub.old %>% tidygraph::as_tibble() %>%
    dplyr::select(c(name, degree))
  n.all <- dplyr::inner_join(n.all, degree.nodes, 
                             by = c("node" = "name"))
  if (nrow(n.all) != old.dim) { stop(call. = TRUE) }
  
  return(list(df.filtered, n.all))
}

##### Function to describe networks
tables.networks <- function(path.to.net, time.point, path.save) {
  #net <- get(load(("../data/intermediate_res_ggm/merged_net.RData")))
  #net <- net$mod_2.2.2.5.5$net
  
  path.data <- "../data/"
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1A")) %>%
    .$feature.data %>%
    dplyr::mutate(CHEBI = as.character(CHEBI)) %>%
    dplyr::mutate(CHEBI = ifelse(grepl("/", CHEBI), 
                                 strsplit(CHEBI, "/")[[1]][1], 
                                 CHEBI)) %>%
    dplyr::select(c(Rvar, var, Class, CHEBI))
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1A")) %>%
    .$feature.data %>%
    dplyr::mutate(CHEBI = as.character(CHEBI)) %>%
    dplyr::mutate(CHEBI = ifelse(grepl("/", CHEBI), 
                                 strsplit(CHEBI, "/")[[1]][1], 
                                 CHEBI)) %>%
    dplyr::select(c(Rvar, var, CHEBI))
  
  net <- get(load(path.to.net))
  if (!is.null(time.point)) {
    net <- net[[paste0("mod_2.2.", time.point, ".5.5")]]$net
  } else {
    net <- net$mod_2.2.2.5.5$net
  }
  
  if (!is.null(time.point)) {
    # Tidy processed but not merged network (i.e., time-specific network)
    net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
      tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.))) %>%
      dplyr::rowwise() %>%
      # Add labels
      dplyr::mutate(label = strsplit(name, "_")[[1]][2])
    
    net <- net %>%
      tidygraph::to_undirected() %>%
      tidygraph::activate(., what = "nodes") %>%
      tidygraph::mutate(label = factor(net.nodes$label)) %>%
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
      dplyr::mutate(name = stringr::str_replace(name, "_.*", "")) %>%
      tidygraph::select(-c(community, degree, type)) %>%
      tidygraph::mutate(group = tidygraph::group_components())
  }
  
  # Add info on chemical class
  all.exposures <- dict.exposure.groups() %>%
    stack() %>% tibble::as_tibble() %>%
    `colnames<-`(c("exposure", "class")) %>%
    dplyr::mutate(class = as.character(class)) %>%
    dplyr::mutate(class = substr(class, 1, nchar(class) - 1)) %>%
    dplyr::mutate(exposure = toupper(exposure))
  net.nodes <- net %>% tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.))) %>%
    dplyr::mutate(label = as.character(label)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(name = ifelse(
      label == "serum metabolome", 
      metabs[metabs$Rvar == name, ]$var, 
      name
    )) %>%
    dplyr::mutate(name = ifelse(
      label == "urinary metabolome", 
      metabu[metabu$Rvar == name, ]$var, 
      name
    )) %>%
    dplyr::mutate(name = stringr::str_remove_all(name, "\"")) %>%
    dplyr::mutate(class = ifelse(
      label == "exposure", all.exposures[all.exposures$exposure == name, ]$class, 
      "none"
    ))
  net <- net %>%
    tidygraph::mutate(name = as.factor(net.nodes$name)) %>%
    tidygraph::mutate(class = as.factor(net.nodes$class)) %>%
    tidygraph::to_undirected() %>%
    #tidygraph::mutate(hub.score = tidygraph::centrality_hub()) %>%
    tidygraph::mutate(degree = as.numeric(tidygraph::centrality_degree()))
  
  # Barplot of `degree` by class
  net.degree <- net %>% tidygraph::as_tibble() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(class = ifelse(
      label == "serum metabolome", 
      metabs[metabs$var == name, ]$Class, 
      ifelse(
        label == "exposure", 
        #all.exposures[all.exposures$exposure == name, ]$class, 
        "exposure", 
        as.character(label))
      )
    )
  net.degree$class <- factor(net.degree$class)
  degree.plt <- net.degree %>%
    dplyr::group_by(class) %>%
    dplyr::top_n(n = 5, wt = degree) %>%
    dplyr::arrange(degree, .by_group = FALSE) %>%
    dplyr::arrange(degree, .by_group = TRUE) %>%
    dplyr::mutate(name = factor(name, 
                                levels = name[order(ave(degree, class, FUN = min), 
                                                    degree)]))
  barplot.degree <- degree.plt %>%
    ggplot2::ggplot() +
    ggplot2::geom_col(aes(x = name, y = degree, 
                          fill = class), 
                      position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(#axis.line = ggplot2::element_blank(), 
                   #axis.title = ggplot2::element_blank(), 
                   #axis.text = ggplot2::element_blank(), 
                   #axis.ticks = ggplot2::element_blank(), 
                   legend.title = element_text(size = 16), 
                   legend.text = element_text(size = 12), 
                   axis.text.x = ggplot2::element_text(angle = 70, hjust = 1)) +
    ggplot2::labs(x = "feature") +
    ggplot2::scale_fill_brewer(palette = "Paired")
  key.save <- ifelse(
    is.null(time.point), "merged", time.point
  )
  ggplot2::ggsave(paste0(path.save, "barplotDegree", "_", key.save, ".png"), 
                  dpi = 720/2, 
                  width = 20, height = 12)
  
  map.label.to.shape <- map.char.to.aes()[[1]]
  map.class.to.col <- map.char.to.aes()[[2]]
  
  ##############################################################################
  # Create table for each chemical class and for full network
  ##############################################################################
  chem.classes <- unique(all.exposures$class)
  all.tables <- list()
  n.all <- list()
  for (chem.class in chem.classes) {
    res.tmp <- .tables.networks(net = net, time.point = time.point, 
                                chem.classes = chem.classes, 
                                chem.class = chem.class)
    all.tables[[chem.class]] <- res.tmp[[1]]
    n.all[[chem.class]] <- res.tmp[[2]]
  }
  full.net <- .tables.networks(net = net, time.point = time.point, 
                               chem.classes = chem.classes, 
                               chem.class = NULL)
  
  # Merge and save table results (description) to file
  ret.nodes <- gtsummary::tbl_summary(n.all %>%
                                        dplyr::bind_rows(.id = "column_label") %>%
                                        dplyr::select(-c(node)), 
                                      by = column_label) %>%
    gtsummary::bold_labels()
  overall <- full.net[[2]] %>% dplyr::select(-c(node)) %>%
    gtsummary::tbl_summary() %>%
    gtsummary::bold_labels()
  overall.nodesAndEdges <- overall %>% gtsummary::add_n()
  ret.nodes <- gtsummary::tbl_merge(list(overall, ret.nodes), 
                                    tab_spanner = FALSE)
  ret <- gtsummary::tbl_summary(all.tables %>%
                                  dplyr::bind_rows(.id = "column_label") %>%
                                  dplyr::mutate(layer.a = factor(layer.a)) %>%
                                  dplyr::mutate(layer.b = factor(layer.b)) %>%
                                  dplyr::mutate(sign = as.numeric(as.character(sign))) %>%
                                  dplyr::select(-c(node.a, node.b, 
                                                   layer.a, layer.b, 
                                                   degree.a, degree.b)), 
                                by = column_label) %>%
    gtsummary::bold_labels()
  overall <- full.net[[1]] %>%
    dplyr::select(-c(node.a, node.b, layer.a, layer.b, 
                     degree.a, degree.b)) %>%
    gtsummary::tbl_summary() %>%
    gtsummary::bold_labels() %>% gtsummary::add_n()
  overall.nodesAndEdges <- gtsummary::tbl_stack(list(overall, 
                                                     overall.nodesAndEdges), 
                                                group_header = c("Edges", "Nodes"))
  ret <- gtsummary::tbl_merge(list(overall, ret), tab_spanner = FALSE)
  combined.res <- gtsummary::tbl_stack(tbls = list(ret, ret.nodes)) %>%
    gtsummary::as_gt() %>%
    gt::gtsave(file = paste0("results/images/desc_net", 
                             ifelse(
                               is.null(time.point), "Merged", time.point
                             ), 
                             "_direct.rtf"))
  
  # Merge and save table results (individual edges) to file
  # Start from full, merged network
  if (is.null(time.point)) {
    ret <- tidy.graph(net) %>%
      dplyr::select(-dplyr::any_of(c("class.x", "class.y", 
                                     "pval.x", "pval.y", "pval", "prob"))) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), round, 3)) %>%
      dplyr::arrange(name.x, name.y, label.x, label.y, 
                     desc(pcor.x), desc(pcor.y)) %>%
      `colnames<-`(c("node.a", "layer.a", "degree.a", 
                     "node.b", "layer.b", "degree.b", 
                     "pcor.a", "qval.a", "sign", "pcor.b", "qval.b")) %>%
      dplyr::select(c(node.a, layer.a, degree.a, 
                      node.b, layer.b, degree.b, 
                      sign, 
                      pcor.a, qval.a, 
                      pcor.b, qval.b))
    colnames(ret) <- c("node.a", "layer.a", "degree.a", 
                       "node.b", "layer.b", "degree.b", 
                       "sign", 
                       "pcor.t1", "qval.t1", "pcor.t2", "qval.t2")
    
    ret <- ret %>%
      dplyr::arrange(node.a, layer.a, node.b, layer.b, 
                     sign, desc(pcor.t1), desc(qval.t1), 
                     desc(pcor.t2), desc(qval.t2)) %>%
      dplyr::select(-c(sign))
    
    # Add genes associated to CpG sites
    ewas.cpgs <- data.table::fread("../data/methylome/results.txt") %>%
      tibble::as_tibble() %>%
      dplyr::select(c("CpG", "Gene"))
    old.dim <- nrow(ret)
    ret <- dplyr::left_join(ret, ewas.cpgs, by = c("node.a" = "CpG"))
    ret <- dplyr::left_join(ret, ewas.cpgs, by = c("node.b" = "CpG"))
    ret <- ret %>% dplyr::distinct()
    if (nrow(ret) != old.dim) { stop(call. = TRUE) }
    if (nrow(dplyr::distinct(ret, node.a, node.b, layer.a, layer.b)) != nrow(ret)) {
      stop(call. = TRUE)
    }
    
    #ret <- ret %>% dplyr::select(-c(column_label))
    ret <- ret %>% dplyr::mutate(id = dplyr::row_number())
    ret[is.na(ret$Gene.x), ]$Gene.x <- ""
    ret[is.na(ret$Gene.y), ]$Gene.y <- ""
    ret <- ret %>%
      tidyr::unite(genes, c(Gene.x, Gene.y), sep = ";")
    ret[ret$genes == ";", ]$genes <- ""
    
    ############################################################################
    # Save sub-table with only mixed interaction, of which one must be an exposure
    ret.mixed <- ret %>%
      dplyr::filter((layer.a == "exposure" & layer.b != "exposure") | 
                      (layer.b == "exposure" & layer.a != "exposure")) %>%
      dplyr::select(-c(id)) %>%
      dplyr::arrange(layer.a, layer.b, node.a, node.b, pcor.t1, pcor.t2)
    pdf(paste0("results/images/mixedEdges_netMerged.pdf"), 
        height = 7, width = 27)
    gridExtra::grid.table(ret.mixed, rows = NULL)
    dev.off()
    write.csv(ret.mixed, file = "results/images/mixedEdges_netMerged.csv", 
              row.names = FALSE)
    ############################################################################
    
    chunk = 20
    n <- nrow(ret)
    r <- rep(1:ceiling(n / chunk), each = chunk)[1:n]
    pdf(paste0("results/images/edges_net", 
               ifelse(
                 is.null(time.point), "Merged", time.point
               ), 
               "_full.pdf"), 
        height = 7, width = 27)
    lapply(split(ret, r), function(d) {
      grid::grid.newpage()
      gridExtra::grid.table(d, rows = NULL)
    })
    dev.off()
  } # End if for table of edges
  
  ##############################################################################
  # Merge and save table results (individual edges) to file
  # Start from full, time-specific networks and filter "top" associations
  ##############################################################################
  if (!is.null(time.point)) {
    ret <- tidy.graph(net) %>%
      dplyr::select(-dplyr::any_of(c("class.x", "class.y", 
                                     "pval", "prob"))) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), round, 3)) %>%
      dplyr::arrange(name.x, name.y, label.x, label.y, 
                     desc(pcor)) %>%
      `colnames<-`(c("node.a", "layer.a", "degree.a", 
                     "node.b", "layer.b", "degree.b", 
                     "pcor", "qval", "sign")) %>%
      dplyr::select(c(node.a, layer.a, degree.a, 
                      node.b, layer.b, degree.b, 
                      sign, 
                      pcor, qval))
    
    ret <- ret %>%
      dplyr::arrange(node.a, layer.a, node.b, layer.b, 
                     sign, desc(pcor), desc(qval)) %>%
      dplyr::select(-c(sign))
    
    # Add genes associated to CpG sites
    ewas.cpgs <- data.table::fread("../data/methylome/results.txt") %>%
      tibble::as_tibble() %>%
      dplyr::select(c("CpG", "Gene"))
    old.dim <- nrow(ret)
    ret <- dplyr::left_join(ret, ewas.cpgs, by = c("node.a" = "CpG"))
    ret <- dplyr::left_join(ret, ewas.cpgs, by = c("node.b" = "CpG"))
    ret <- ret %>% dplyr::distinct()
    if (nrow(ret) != old.dim) { stop(call. = TRUE) }
    if (nrow(dplyr::distinct(ret, node.a, node.b, layer.a, layer.b)) != nrow(ret)) {
      stop(call. = TRUE)
    }
    
    ret <- ret %>% dplyr::mutate(id = dplyr::row_number())
    ret[is.na(ret$Gene.x), ]$Gene.x <- ""
    ret[is.na(ret$Gene.y), ]$Gene.y <- ""
    ret <- ret %>%
      tidyr::unite(genes, c(Gene.x, Gene.y), sep = ";")
    ret[ret$genes == ";", ]$genes <- ""
    
    ############################################################################
    # Save sub-table with only mixed interaction, of which one must be an exposure
    ############################################################################
    ret.mixed <- ret %>%
      dplyr::filter((layer.a == "exposure" & layer.b != "exposure") | 
                      (layer.b == "exposure" & layer.a != "exposure")) %>%
      dplyr::select(-c(id, genes)) %>%
      dplyr::arrange(desc(abs(pcor)), layer.a, layer.b, node.a, node.b)
    write.csv(ret.mixed, file = paste0("results/images/mixedEdges_net", 
                                       time.point, ".csv"), 
              row.names = FALSE)
    
    chunk = 20
    n <- nrow(ret.mixed)
    r <- rep(1:ceiling(n / chunk), each = chunk)[1:n]
    pdf(paste0("results/images/edges_net", 
               ifelse(
                 is.null(time.point), "Merged", time.point
               ), 
               "_mixed.pdf"), 
        height = 7, width = 27)
    lapply(split(ret.mixed, r), function(d) {
      grid::grid.newpage()
      gridExtra::grid.table(d, rows = NULL)
    })
    dev.off()
  } # End if for table of edges for time-specific networks
  
  return(overall.nodesAndEdges)
}
