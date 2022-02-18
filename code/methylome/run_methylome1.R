setwd("/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/")
source("code/methylome/methylome.R")

meth <- load.methylome(filter.time = list(
  perform = FALSE, time.point = NULL
))
res <- methylome.ewaff(meth = meth, time.point = 1, perform.adj = TRUE)
