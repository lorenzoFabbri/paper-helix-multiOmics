library(WGCNA)

source("code/methylome/methylome.R")

meth <- load.methylome(filter.time = list(
  perform = FALSE, time.point = NULL
))

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 18)

datExpr <- meth %>%
  dplyr::filter(grepl("_1A", SampleID)) %>%
  dplyr::select(-c(SampleID, HelixID)) %>%
  as.matrix()

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
print(gsg$allOK)
sampleTree = hclust(dist(datExpr), method = "average")
gc()

bwnet = WGCNA::blockwiseModules(datExpr, 
                                maxBlockSize = 1000, minModuleSize = 30, 
                                power = 6, TOMType = "unsigned", 
                                reassignThreshold = 0, mergeCutHeight = 0.25, 
                                numericLabels = TRUE, 
                                saveTOMs = FALSE, saveTOMFileBase = "methyl", 
                                verbose = 7, nThreads = 18)
