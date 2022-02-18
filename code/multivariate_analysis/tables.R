# Author: Lorenzo Fabbri
# Script containing functions to create "publication-ready" tables with results

library(readr)
library(tibble)
library(xtable)
library(gridExtra)
library(ggplot2)
library(knitr)
library(png)
library(grid)
library(stringr)
library(purrr)


# Function to make report w/ important plots arranged by time point
make.report.by.time <- function(path.all.plots, pcs) {
  len.code <- 6 # model.code
  name.plots.compIcompJ <- paste0("comp", 1:(pcs - 1), 
                                  "comp", 2:pcs, 
                                  "_[X,Y]_mod_.*\\.png")
  name.plots.loadings <- paste0("loadings_comp", 1:pcs, "_mod_")
  name.plots.r2 <- paste0("r2_summary_comp", 1:pcs, "_mod_")
  name.plots.circos <- paste0("corr_comp", 1:(pcs - 1), 
                              "comp", 2:pcs, 
                              "_mod_")
  name.plots <- c(name.plots.compIcompJ, name.plots.loadings, 
                  "q2_summary_mod_", "q2_total_mod_", 
                  name.plots.r2, name.plots.circos)
  
  # Format mode.code: scaling.type . adjust.type . model.type . time.point . omic.type . exp.type
  # Iterate over type of plot
  for (name.plot in name.plots) {
    if (name.plot %in% name.plots.compIcompJ) {
      # This is just to that we have X and Y in the same file
      file.names <- list.files(path.all.plots, 
                               pattern = name.plot, 
                               full.names = TRUE)
      pc_i <- strsplit(name.plot, "comp")[[1]][2]
      pc_j <- strsplit(name.plot, "comp")[[1]][3] %>% strsplit(., "_") %>%
        unlist() %>% .[1]
      name.plot <- paste0("comp", pc_i, 
                          "comp", pc_j, 
                          "_XY_mod_")
    } else {
      file.names <- list.files(path.all.plots, 
                               pattern = glob2rx(paste0(name.plot, "*.png")), 
                               full.names = TRUE)
    }
    
    # Extract model.code for specific type of plot
    model.codes <- list()
    for (path in file.names) {
      tmp <- substr(path, nchar(path) - (len.code + 4) + 1, nchar(path)) %>%
        strsplit(., split = "[.]") %>% unlist() %>% .[[1]] %>%
        strsplit(., split = "") %>% unlist()
      model.codes <- append(model.codes, list(tmp))
    }
    model.codes <- as.data.frame(matrix(unlist(model.codes), 
                                        ncol = len.code, byrow = TRUE))
    model.codes <- model.codes %>% dplyr::mutate_all(as.factor)
    if (len.code == 4) {
      colnames(model.codes) <- c("scaling.type", "model.type", "time.point", "omic.type")
      to.split.by <- list("scaling.type", "model.type", "omic.type")
      label.plot <- c("s", "m", "o")
    } else if (len.code == 5) {
      colnames(model.codes) <- c("scaling.type", "adj.type", "model.type", 
                                 "time.point", "omic.type")
      to.split.by <- list("scaling.type", "adj.type", "model.type", "omic.type")
      label.plot <- c("s", "a", "m", "o")
    } else if (len.code == 6) {
      colnames(model.codes) <- c("scaling.type", "adj.type", "model.type", 
                                 "time.point", "omic.type", "exp.type")
      to.split.by <- list("scaling.type", "adj.type", "model.type", 
                          "omic.type", "exp.type")
      label.plot <- c("s", "a", "m", "o", "e")
    } else {
      stop(call. = TRUE)
    }
    
    # Grids by time point: separate files for scaling.type, adj.type, 
    # model.type, omic.type, exp.type
    model.codes.group <- model.codes %>%
      dplyr::group_by_at(setdiff(names(.), "time.point")) %>%
      dplyr::group_split(.)
    
    for (grp in model.codes.group) {
      idxs <- grp %>% tidyr::unite(idxs, names(grp), sep = "") %>%
        dplyr::pull() %>% as.vector()
      
      idxs.collapsed <- grep(paste0(idxs, collapse = "|"), 
                             file.names, value = TRUE)
      plots <- lapply(idxs.collapsed, function(x) {
        imgs <- as.raster(png::readPNG(x))
        grid::rasterGrob(imgs, interpolate = FALSE, 
                         width = unit(10, "in"), height = unit(6, "in"))
      })
      
      label.plot.all <- strsplit(idxs[1], "")[[1]][-2] %>%
        paste0(., label.plot, collapse = "")
      ggplot2::ggsave(filename = paste0(name.plot, label.plot.all, 
                                        "_", "all.jpg"),
                      plot = gridExtra::marrangeGrob(grobs = plots,
                                                     ncol = 2,
                                                     nrow = ceiling(length(plots) / 2)),
                      device = "jpg",
                      path = path.all.plots, dpi = 270, 
                      width = unit(20, "in"), height = unit(12, "in"))
    } # End loop over grouped dataframe
    
  } # End loop over name.plots
  
}

# Function to read file containing unstructured results
read.file.res <- function(file.path) {
  res <- readr::read_lines(file = file.path)
  
  # Create df from lines
  model.code   <- res[1]
  omic.type    <- res[2]
  time.point   <- res[3]
  model.type   <- ifelse(res[4] == "FALSE", "PLS", "sPLS")
  adjusted     <- res[5]
  type.scaling <- res[6]
  group.scaling <- res[7]
  max.iter     <- res[8]
  nrepeat      <- res[9]
  nfolds       <- res[10]
  exp.dim      <- res[11]
  omics.dim    <- res[12]
  if (res[4] == FALSE) { # PLS
    prop.var.expl.X <- res[13]
    prop.var.expl.Y <- res[14]
    mean.Q2.total   <- res[17]
    sd.Q2.total     <- res[18]
    min.Q2.total    <- res[19]
    max.Q2.total    <- res[20]
  } else { # sPLS
    ncomp <- res[13]
    keepX <- res[14]
    keepY <- res[15]
    prop.var.expl.X <- res[16]
    prop.var.expl.Y <- res[17]
    mean.Q2.total   <- res[20]
    sd.Q2.total     <- res[21]
    min.Q2.total    <- res[22]
    max.Q2.total    <- res[23]
  }
  
  prop.var.expl.X <- lapply(prop.var.expl.X, 
                            function(x) {str2lang(x) %>% eval(.) %>%
                                as.vector(.) %>% round(., digits = 3) %>%
                                unlist(.)})
  prop.var.expl.Y <- lapply(prop.var.expl.Y, 
                            function(x) {str2lang(x) %>% eval(.) %>%
                                as.vector(.) %>% round(., digits = 3) %>%
                                unlist(.)})
  mean.Q2.total <- lapply(mean.Q2.total, 
                          function(x) {str2lang(x) %>% eval(.) %>%
                              as.vector(.) %>% round(., digits = 3) %>%
                              unlist(.)})
  sd.Q2.total <- lapply(sd.Q2.total, 
                        function(x) {str2lang(x) %>% eval(.) %>%
                            as.vector(.) %>% round(., digits = 3) %>%
                            unlist(.)})
  min.Q2.total <- lapply(min.Q2.total, 
                         function(x) {str2lang(x) %>% eval(.) %>%
                             as.vector(.) %>% round(., digits = 3) %>%
                             unlist(.)})
  max.Q2.total <- lapply(max.Q2.total, 
                         function(x) {str2lang(x) %>% eval(.) %>%
                             as.vector(.) %>% round(., digits = 3) %>%
                             unlist(.)})
  
  exp.type <- strsplit(model.code, "_")[[1]][2] %>% strsplit(., "") %>%
    unlist() %>% .[6] %>% as.character()
  exposure.type <- list("1" = "phenols", 
                        "2" = "phthalates.low", "3" = "phthalates.high", 
                        "4" = "pesticides")
  
  ret <- tibble::tibble(model.code, omic.type, exposure.type[exp.type], 
                        time.point, 
                        model.type, adjusted, 
                        type.scaling, group.scaling, 
                        #max.iter, nrepeat, nfolds, 
                        exp.dim, omics.dim)
  if (res[4] == FALSE) {
    ret <- ret %>% tibble::add_column(prop.var.expl.X, 
                                      prop.var.expl.Y, 
                                      mean.Q2.total, sd.Q2.total)
                                      #min.Q2.total, max.Q2.total)
  } else {
    ret <- ret %>% tibble::add_column(ncomp, keepX, keepY, 
                                      prop.var.expl.X, 
                                      prop.var.expl.Y, 
                                      mean.Q2.total, sd.Q2.total, 
                                      min.Q2.total, max.Q2.total)
  }
  
  return(ret)
}

generate.table <- function(files.path) {
  # Iterate over all files in directory and create df
  file.names <- list.files(files.path, pattern = "*.txt", 
                           full.names = TRUE)
  dfs <- lapply(file.names, read.file.res)
  dfs.merged <- dplyr::bind_rows(dfs)
  options(knitr.kable.NA = "")
  table.res <- knitr::kable(dfs.merged, format = "pipe", digits = 3)
  
  # Add all tables to PDF document
  # for (df in dfs) {
  #   tg <- gridExtra::tableGrob(df)
  #   h  <- grid::convertHeight(sum(tg$heights), "mm", TRUE)
  #   w  <- grid::convertWidth(sum(tg$widths), "mm", TRUE)
  #   
  #   write.path.tmp <- paste(write.path, df$model.code, ".pdf")
  #   ggplot2::ggsave(write.path.tmp, tg, 
  #                   width = w, height = h, units = 'mm')
  # }
  
  return(dfs.merged)
}
