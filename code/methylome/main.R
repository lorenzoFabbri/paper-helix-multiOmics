# Call to functions to load the results of the EWASs (by time point)
# and process them to find the common CpG sites for the GGMs

source("code/methylome/methylome.R")

res.ewaff <- load.res.ewaff(path.res = "results/ewaff_win_sva/")
processed.ewaff <- process.res.ewaff(df1 = res.ewaff[[1]], 
                                     df2 = res.ewaff[[2]], 
                                     key.save = "sva", 
                                     threshold.fdr = 0.2)
