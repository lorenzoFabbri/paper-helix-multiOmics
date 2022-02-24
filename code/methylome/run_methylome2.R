setwd("/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/")
source("code/methylome/methylome.R")

options(mc.cores = 25)

meth <- load.methylome(filter.time = list(
  perform = FALSE, time.point = NULL
))
res <- methylome.ewaff(meth = meth, time.point = 2, perform.adj = TRUE)
