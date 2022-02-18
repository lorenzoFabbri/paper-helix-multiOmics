# Author: Lorenzo Fabbri

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(broom)
library(stats)
library(ComplexHeatmap)
library(labelled)
library(igraph)
library(gsubfn)
library(stringr)
library(ggdendro)
library(dendextend)
library(readxl)
library(MetaboAnalystR)
library(ggplot2)

source("./code/multivariate_analysis.R")

#################### PART I ####################
# Regress each exposure against each PC of -omics
regress.exps.pcOmes <- function(exposome, omics, metadata, 
                                number.pcs, scale.perf) {
  # Pre-process datasets
  tmp      <- preprocess.data(exposome, omics, metadata, 
                              scale.perf = scale.perf)
  exposome <- tmp[[1]]
  omic     <- tmp[[2]]
  metadata <- tmp[[3]]
  
  # Perform PCA on -omics and extract PCs
  pca.res <- FactoMineR::PCA(omic, ncp = number.pcs, 
                             scale. = FALSE, graph = FALSE)
  principal.components <- as.data.frame(pca.res$ind$coord)
  
  # Loop to regress each exposure against each PC of -omics
  mods.res <- list()
  for (exposure in colnames(exposome)) {
    for (p.comp in colnames(principal.components)) {
      dat <- as.data.frame(cbind(as.matrix(exposome[, exposure]), 
                                 as.matrix(principal.components[, p.comp])))
      colnames(dat) <- c(exposure, p.comp)
      
      mod.res <- lm(as.formula(paste(exposure, " ~ ", p.comp)), 
                    data = dat)
      mods.res <- append(mods.res, list(mod.res))
    }
  }
  
  # Adjust models
  ret <- data.frame()
  i <- 1
  col.names.exp <- rep(colnames(exposome), each = number.pcs)
  for (mod in mods.res) {
    mod.tmp <- broom::tidy(mod)[2, ] %>%
      dplyr::mutate(exposure = col.names.exp[i]) %>%
      dplyr::mutate(n_obs = stats::nobs(mod))
    
    i <- i + 1
    ret <- rbind(ret, mod.tmp)
  }
  
  ret.adj <- ret %>%
    group_by(term) %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr"))
  
  return(ret.adj)
}

#################### PART II ####################
# Hierarchical clustering of Spearman correlation for exposures
hclust.correlation.exposome <- function(exposome, omics, metadata, 
                                        scale.perf = scale.perf, threshold) {
  
  # Family group for exposures
  labs <- exposome %>%
    dplyr::select(-c(HelixID, SubjectID)) %>%
    dplyr::select(where(is.numeric)) %>%
    labelled::var_label(.) %>%
    unlist() %>% as.data.frame() %>% rownames_to_column(., "exposure")
  colnames(labs) <- c("exposure", "family")
  
  tmp      <- preprocess.data(exposome, omics, metadata, 
                              scale.perf = scale.perf)
  exposome <- tmp[[1]]
  
  # For the sake of visualization, rename the exposures
  # Fucking stupid to use a for loop...
  old.colnames <- colnames(exposome)
  new.colnames <- list()
  for (el in old.colnames) {
    tmp.el <- str_split(el, "\\.")[[1]][2]
    tmp.el <- str_split(tmp.el, "\\_")[[1]][2]
    new.colnames <- append(new.colnames, tmp.el)
  }
  #colnames(exposome) <- new.colnames
  
  # Compute Spearman rank correlation coefficients
  spearman <- cor(exposome, method = "spearman")
  
  plt <- ComplexHeatmap::Heatmap(spearman, 
                                 cluster_rows = T, #left_annotation = ha, 
                                 show_row_dend = T, show_column_dend = F, 
                                 show_column_names = T, show_row_names = T, 
                                 use_raster = T, 
                                 name = 'Heatmap Spearman correlation', 
                                 heatmap_legend_param = list(at = (c(-1, 0 , 1))))
  print(plt)
  
  # Convert data into correlation network
  col <- c("gray50", "tomato", "gold")
  
  corr.graph <- as.data.frame(spearman)
  corr.graph <- as.matrix(corr.graph)
  corr.graph <- igraph::graph.adjacency(corr.graph, 
                                        weighted = TRUE, mode = "undirected", diag = FALSE)
  corr.graph <- igraph::delete.edges(corr.graph, igraph::E(corr.graph)[ abs(weight) < threshold ])
  
  igraph::V(corr.graph)$short  <- new.colnames
  igraph::V(corr.graph)$family <- as.factor(labs$family)
  igraph::V(corr.graph)$color  <- col[as.integer(igraph::V(corr.graph)$family)]
  plot(corr.graph, 
       vertex.size = 27, vertex.label = unlist(igraph::V(corr.graph)$short), 
       edge.color = ifelse(spearman > 0, "blue", "red"), 
       edge.width = abs(igraph::E(corr.graph)$weight) * 8)
  legend("topleft", 
         legend = levels(as.factor(labs$family)), fill = col)
  
  # Community detection
  communities.corr <- igraph::cluster_edge_betweenness(corr.graph, 
                                                       directed = FALSE)
  plot(communities.corr, corr.graph, 
       col = igraph::V(corr.graph)$color, 
       vertex.label = unlist(igraph::V(corr.graph)$short))
  legend("topleft", 
         legend = levels(as.factor(labs$family)), fill = col)
}

# Clustering of -omics data
hclust.omics <- function(exposome, omics, metadata, 
                         scale.perf = scale.perf, threshold) {
  # Data pre-processing
  tmp      <- preprocess.data(exposome, omics, metadata, 
                              scale.perf = scale.perf)
  omics <- tmp[[2]]
  
  # Clustering
  d    <- dist(t(omics), method = "euclidean")
  hc   <- hclust(d,   method = "ward.D")
  dend <- as.dendrogram(hc)
  labels(dend) <- colnames(omics[order.dendrogram(dend)])
  
  # Cut tree
  leaves <- threshold
  dendextend::cutree(dend, k = leaves)[order.dendrogram(dend)]
  dend.col <- color_branches(dend,   k = leaves)
  dend.col <- color_labels(dend.col, k = leaves)
  print(plot(dend.col))
  
  # Pathway enrichment analysis
  perform.pathEnrichment.omics(dend, omics, threshold)
  
  return(dend)
}

perform.pathEnrichment.omics <- function(dend, omics, threshold) {
  PATH = '/run/user/1000/gvfs/smb-share:server=fs01.isglobal.lan,share=data/ws_helix/HELIX_preproc/metabolome/Final_data/'
  type.omics <- ifelse(grepl("blood", deparse(substitute(omics))), "serum", "urine")
  annotations <- paste0("Metabolome_", type.omics, "_annotation_v2.xlsx")
  file.annotation <- paste0(PATH, annotations)
  
  # Cut tree
  communities <- dendextend::cutree(dend, k = threshold)[order.dendrogram(dend)]
  communities.df <- as.data.frame(unlist(communities)) %>%
    rownames_to_column()
  colnames(communities.df) <- c("metabolite", "community")
  
  # Add CHEBI and KEGG codes to metabolites
  annotations <- readxl::read_xlsx(file.annotation, sheet = 1)
  
  communities.codes <- communities.df %>%
    dplyr::inner_join(., annotations[, c("Rvar", "CHEBI", "KEGG")], 
                      by = c("metabolite" = "Rvar")) %>%
    dplyr::filter(CHEBI != 'NA' & KEGG != 'NA') %>%
    dplyr::filter(!(CHEBI %>% str_detect("/"))) %>%
    dplyr::filter(!(KEGG %>% str_detect("/")))
  
  ## Use MetaboAnalyst to perform pathway enrichment analysis (PEA)
  # Internal function to perform PEA
  perform.pea <- function(communities, type, cluster) {
    mSet <- MetaboAnalystR::InitDataObjects("conc", "pathora", FALSE)
    mSet <- MetaboAnalystR::Setup.MapData(mSet, communities$KEGG)
    mSet <- MetaboAnalystR::CrossReferencing(mSet, "hmdb_kegg")
    mSet <- MetaboAnalystR::CreateMappingResultTable(mSet)
    
    # Detailed search of missing metabolites (i.e., no match)
    not.matched <- as.data.frame(mSet$dataSet$map.table) %>%
      filter(Match == 'NA')
    for (q in not.matched$Query) {
      q.metabolite <- communities[communities$KEGG == q, ]$metabolite
      if (startsWith(q.metabolite, c("X"))) {
        q.metabolite <- substring(q.metabolite, 2)
      }
      
      mSet <- MetaboAnalystR::PerformDetailMatch(mSet, q.metabolite)
      mSet <- MetaboAnalystR::GetCandidateList(mSet)
    }
    
    # Ora
    mSet <- MetaboAnalystR::SetKEGG.PathLib(mSet, "hsa", "current")
    mSet <- MetaboAnalystR::SetMetabolomeFilter(mSet, F)
    mSet <- MetaboAnalystR::CalculateOraScore(mSet, "rbc", "hyperg")
    
    # Plot results
    mSet.results <- as.data.frame(mSet$api$ora.results)
    mSet.plot <- ggplot(mSet.results, aes(x = Impact, y = `-log(p)`)) +
      geom_point(aes(color = `-log(p)`, size = Impact)) +
      scale_fill_distiller(palette = "Blues", direction = -1) +
      geom_text(label = ifelse(mSet.results$`-log(p)` > 2 & mSet.results$Impact > 0.1, 
                               rownames(mSet.results), "")) +
      ggtitle(paste(type, ": ", cluster))
    print(mSet.plot)
  }
  
  # Loop over found communities and perform PEA
  for (community in unique(communities.codes$community)) {
    perform.pea(communities = communities.codes[communities.codes$community == community, ], 
                type = deparse(substitute(omics)), cluster = as.character(community))
  }
}
