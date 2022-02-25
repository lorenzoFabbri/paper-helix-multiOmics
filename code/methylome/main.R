source("code/methylome/methylome.R")

res.ewaff <- load.res.ewaff(path.res = "results/ewaff_win_sva/")
processed.ewaff <- process.res.ewaff(df1 = res.ewaff[[1]], 
                                     df2 = res.ewaff[[2]], 
                                     key.save = "sva", 
                                     threshold.fdr = 0.2)

res.ewaff <- load.res.ewaff(path.res = "results/ewaff_win/")
processed.ewaff <- process.res.ewaff(df1 = res.ewaff[[1]], 
                                     df2 = res.ewaff[[2]], 
                                     key.save = "tmp", 
                                     threshold.fdr = 0.2)
