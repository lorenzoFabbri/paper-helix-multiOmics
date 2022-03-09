# "Script" that contains all the dictionaries and mappings used
# Author: Lorenzo Fabbri

library(webchem)

# Returns information for MeSH analyses
mesh.analyses <- function() {
  # Files for MeSH analysis
  mesh.res <- get(load(url("https://github.com/barupal/ChemRICH/blob/master/cidmesh_smiles_fpbit.RData?raw=true"))) %>%
    tibble::as_tibble()
  tree.res <- get(load(url("https://github.com/barupal/ChemRICH/blob/master/treenames.df.RData?raw=true"))) %>%
    tibble::as_tibble()
  
  # Omic variables
  params <- list(time.point = 1)
  path.data <- "../data/"
  prot <- readRDS(file = paste0(path.data, "proteome_1", 
                                ifelse(params$time.point == 1, 
                                       "A", "B"))) %>% .$feature.data
  metabs <- readRDS(file = paste0(path.data, "metabSerum_1", 
                                  ifelse(params$time.point == 1, 
                                         "A", "B"))) %>% .$feature.data
  metabu <- readRDS(file = paste0(path.data, "metabUrine_1", 
                                  ifelse(params$time.point == 1, 
                                         "A", "B"))) %>% .$feature.data
  
  # Convert compounds' codes to PubChem CID
  metabs <- metabs %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pcid = webchem::cts_convert(query = paste0("CHEBI:", CHEBI), 
                                              from = "chebi", to = "pubchem cid") %>%
                    unname() %>% .[[1]] %>% .[1])
  metabu <- metabu %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pcid = webchem::cts_convert(query = paste0("CHEBI:", CHEBI), 
                                              from = "chebi", to = "pubchem cid") %>%
                    unname() %>% .[[1]] %>% .[1])
  
  # Add information from DBs to our datasets
  metabs <- metabs %>%
    dplyr::left_join(mesh.res, by = c("pcid" = "CID")) %>%
    dplyr::left_join(tree.res, by = c("MESSTREE" = "MESHTREE")) %>%
    dplyr::select(-c(MeSHUID, CIDCount, PMIDCount, fpbit, Rvar.log, 
                     ULOQ, LLOQ, Indicators_extra))
  metabu <- metabu %>%
    dplyr::left_join(mesh.res, by = c("pcid" = "CID")) %>%
    dplyr::left_join(tree.res, by = c("MESSTREE" = "MESHTREE")) %>%
    dplyr::select(-c(MeSHUID, CIDCount, PMIDCount, fpbit, Rvar.log))
}

# Returns df (from tidygraph's graph) with labels instead of indices in tibble of edges
tidy.graph <- function(net) {
  node.attributes <- net %>% tidygraph::activate(., what = "nodes") %>%
    tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
  edge.attributes <- net %>% tidygraph::activate(., what = "edges") %>%
    tidygraph::as_tibble()
  
  df <- dplyr::inner_join(node.attributes, edge.attributes, 
                          by = c("idx" = "from"))
  df <- dplyr::inner_join(node.attributes, df, 
                          by = c("idx" = "to"))
  df <- df %>%
    dplyr::select(-dplyr::any_of(c("idx", "idx.y", "type.x", "type.y", 
                                   "group.x", "group.y")))
  
  if (dim(df)[1] != dim(edge.attributes)[1]) { stop(call. = TRUE) }
  
  return(df)
}

# Returns named vectors for mapping classes and labels to colors and shapes
map.char.to.aes <- function() {
  map.label.to.shape <- c(
    exposure = 23, 
    methylome = 21, 
    `serum metabolome` = 25, 
    `urinary metabolome` = 24, 
    proteome = 22
  )
  map.class.to.col <- c(
    phthalates.high = "#9E0142", 
    phthalates.low = "#D53E4F", 
    pesticides = "#F46D43", 
    phenols = "#FDAE61", 
    glycerophospholipids = "#FEE08B", 
    none = "#5E4FA2", 
    aminoacids = "#E6F598", 
    acylcarnitines = "#ABDDA4", 
    sphingolipids = "#66C2A5", 
    biogenicamines = "#3288BD", 
    sugars = "#FFFFBF"
  )
  map.layer.to.col <- c(
    exposure = "#CCCCCC", 
    methylome = "#66CCFF", 
    `serum metabolome` = "#FF3333", 
    `urinary metabolome` = "#FFCC00", 
    proteome = "#FFCCCC"
  )
  
  return(list(map.label.to.shape, map.class.to.col, map.layer.to.col))
}

# Returns lists of variables for each -omic layer
vars.detailed <- function() {
  load("./data/env")
  
  # Exposures
  exposures <- base::union(colnames(exp_metab_blood), colnames(exp_metab_blood2))
  exposures <- exposures[!(exposures %in% c("SampleID", "HelixID"))] %>%
    as.data.frame() %>% `colnames<-`("name") %>%
    dplyr::rowwise() %>%
    dplyr::transmute(name = ifelse(grepl("log.hcp_", name), 
                                   strsplit(name, split = "_")[[1]][2], 
                                   name)) %>%
    dplyr::mutate(type = "exp")
  
  # Omics
  serum.met <- union(colnames(omic_metab_blood), colnames(omic_metab_blood2))
  serum.met <- serum.met[!(serum.met %in% c("SampleID", "HelixID"))] %>%
    as.data.frame() %>% `colnames<-`("name") %>%
    dplyr::rowwise() %>%
    dplyr::transmute(name = ifelse(grepl("log.", name), 
                                   strsplit(name, split = "log.")[[1]][2] %>%
                                     paste0(., "_metab_blood"), 
                                   name %>%
                                     paste0(., "_metab_blood"))) %>%
    dplyr::mutate(type = "serum")
  urinary.met <- union(colnames(omic_metabun), colnames(omic_metabun2))
  urinary.met <- urinary.met[!(urinary.met %in% c("SampleID", "HelixID"))] %>%
    as.data.frame() %>% `colnames<-`("name") %>%
    dplyr::rowwise() %>%
    dplyr::transmute(name = ifelse(grepl("log.", name), 
                                   strsplit(name, split = "log.")[[1]][2] %>%
                                     paste0(., "_metabun"), 
                                   name %>%
                                     paste0(., "_metabun"))) %>%
    dplyr::mutate(type = "urine")
  proteome <- union(colnames(omic_proteome), colnames(omic_proteome2))
  proteome <- proteome[!(proteome %in% c("SampleID", "HelixID"))] %>%
    as.data.frame() %>% `colnames<-`("name") %>%
    dplyr::rowwise() %>%
    dplyr::transmute(name = ifelse(grepl("log.", name), 
                                   strsplit(name, split = "log.")[[1]][2] %>%
                                     paste0(., "_proteome"), 
                                   name %>%
                                     paste0(., "_proteome"))) %>%
    dplyr::mutate(type = "prot")
  
  ret <- bind_rows(exposures, serum.met, urinary.met, proteome)
  return(ret)
}

# Returns list of covariates for adjustments
list.covariates <- function() {
  covariates.char <- c("e3_sex.x", "cohort.x", "season")
  covariates.num <- c("age_sample_years.x", "zBMI", 
                      "hs_dift_mealblood_imp")
  covariates <- list(covariates.char = covariates.char, 
                     covariates.num = covariates.num)
  
  return(covariates)
}

# Mapping from parameters (strings) to numerical code for model.code
params.to.model.code <- function(model.type) {
  code.scaling.type <- list(range = 1, autoscale = 2, none = 3)
  code.adj.type <- list(variables = 1, residual = 2, none = 4)
  
  if (model.type == "xPLS") {
    code.model.type <- list(pls = 1, spls = 2, sgpls = 3)
  } else if (model.type == "GMM") {
    code.model.type = list(pearson = 1, spearman = 2, ppcor = 3)
  } else {
    print("Method not supported. Options are `xPLS` or `GMM`.")
  }
  #code.time.type
  code.omic.type <- c(metabun = 1, metab_blood = 2, proteome = 3, 
                      methylome = 4, all = 5)
  code.exposure.type <- list(phenols = 1, 
                             phthalates.low = 2, phthalates.high = 3, 
                             pesticides = 4, all = 5)
  
  ll <- list(
    code.scaling.type = code.scaling.type, 
    code.adj.type = code.adj.type, 
    code.model.type = code.model.type, 
    code.omic.type = code.omic.type, 
    code.exposure.type = code.exposure.type
    )
  return(ll)
}

# From `SampleID` to `HelixID`
from.sample.to.helix <- function(sample.id.col) {
  sapply(sample.id.col, function(x) {
    substr(x, 1, nchar(x) - 3) %>% unname() %>% unlist()
  })
}

# From code of exposure to actual name
from.code.to.exposure <- function(code.exp) {
  # Tidy code exposure to retrieve only actual exposure name
  if (startsWith(code.exp, "log.")) {
    code.exp <- strsplit(code.exp, "log.hcp_,|_") %>%
      unlist() %>% .[2]
  }
  
  dict.exposure.names <- list(
    # Phenols (7)
    mepa = list(name = "methyl-paraben", id = "7456"), 
    etpa = list(name = "ethyl-paraben", id = "8434"), 
    prpa = list(name = "propyl-paraben", id = "7175"), 
    bupa = list(name = "n-butyl-paraben", id = "7184"), 
    bpa = list(name = "bpa", id = "6623"), 
    trcs = list(name = "triclosan", id = "5564"), 
    oxbe = list(name = "oxybenzone", id = "4632"), 
    # Phthalates (10)
    mep = list(name = "monoethyl-phthalate", id = "75318"), 
    mibp = list(name = "mono-iso-butyl-phthalate", id = "92272"), 
    mnbp = list(name = "mono-n-butyl-phthalate", id = "8575"), 
    mehp = list(name = "mono-2-ethylhexyl-phthalate", id = "20393"), 
    mehhp = list(name = "mono-2-ethyl-5-hydroxyhexyl-phthalate", id = "170295"), 
    meohp = list(name = "mono-2-ethyl-5-oxohexyl-phthalate", id = "119096"), 
    mecpp = list(name = "mono-2-ethyl-5-carboxypentyl-phthalate", id = "148386"), 
    mbzp = list(name = "mono-benzyl-phthalate", id = "31736"), 
    ohminp = list(name = "mono-4-methyl-7-hydroxyoctyl-phthalate", id = NA), 
    oxominp = list(name = "mono-4-methyl-7-oxooctyl-phthalate", id = "131698778"), 
    # DAPs (should be 6)
    dmtp = list(name = "dimethyl-thiophosphate", id = "168140"), 
    dep = list(name = "diethyl-phosphate", id = "6781"), 
    detp = list(name = "diethyl-thiophosphate", id = "3683036")
  )
  
  if (code.exp == "all") {
    return(dict.exposure.names)
  }
  
  return(dict.exposure.names[code.exp])
}

# Dictionary of exposure groups/families
dict.exposure.groups <- function() {
  exposures <- c(phenols = c("mepa", "etpa", "prpa", "bupa", 
                                "bpa", "trcs", "oxbe"), 
                    phthalates.low = c("mep", "mibp", "mnbp"), 
                    phthalates.high = c("mehp", "mehhp", "meohp", "mecpp", "mbzp", 
                                        "ohminp", "oxominp"), 
                    pesticides = c("dmp", "dmdtp", "dmtp", "dep", "detp"))
  return(exposures)
}
