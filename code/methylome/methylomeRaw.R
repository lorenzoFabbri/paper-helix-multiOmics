library(tidyverse)
library(minfi)

load("../data/methylome/methylome_panel_ComBatSlide_6cells_v4.Rdata")
dat <- get("methylome_panel_ComBatSlide_6cells")
#load("../data/methylome/methylome_panel_residuals2beta_v4.Rdata")
#dat <- get("methylome_panel_residuals2beta")

minfi::densityPlot(dat@assays@.xData$data$Beta)
print(minfi::preprocessMethod(dat))

# Processing
dat <- minfi::dropLociWithSnps(dat, snps = c("SBE", "CpG"), maf = 0)
gc()

# Add metadata
SampleID <- dat@colData$SampleID
HelixID <- dat@colData$HelixID
gc()

dat <- dat@assays@.xData$data$Beta %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(., var = "SampleID")
gc()
dat$HelixID <- HelixID

data.table::fwrite(dat, 
                   file = "../data/methylome/methylome_panel_ComBatSlide_6cells_v4.csv", 
                   nThread = 6)
