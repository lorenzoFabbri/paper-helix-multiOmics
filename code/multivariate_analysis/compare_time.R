# Author: Lorenzo Fabbri
# Script containing functions to compare results by time point

library(gridExtra)

# Function to plot compX of t1 vs. compX of t2 for all components 
# for all models stored in given directory
plot.comparison.variates.wrap <- function(path.models) {
  idx.time.point <- 4
  
  # Iterate over all directories containing saved images
  file.names <- list.files(path.models, 
                           full.names = TRUE)
  models1 <- list()
  models2 <- list()
  for (dir in file.names) {
    code <- strsplit(dir, "/")[[1]][6] %>% strsplit(., "_") %>% unlist(.) %>%
      .[2] %>% strsplit(., "") %>% unlist()
    
    if (code[idx.time.point] == "1") {
      models1 <- append(models1, list(dir))
    } else if (code[idx.time.point] == "2") {
      models2 <- append(models2, dir)
    } else {
      stop(call. = TRUE)
    }
  } # End loop over list directories
  
  for (idx in 1:length(models1)) {
    model1 <- get(load(paste0(models1[idx], "/ret.Rdata")))$model
    model2 <- get(load(paste0(models2[idx], "/ret.Rdata")))$model
    params <- get(load(paste0(models1[idx], "/pls.params.Rdata")))
    
    # Run function to plot results
    plot.comparison.variates(model1, model2, params);
  }
  
} # End function

# Function to plot compX of t1 vs. compX of t2 for all components
plot.comparison.variates <- function(model1, model2, params1) {
  ncomp <- model1$ncomp
  model.code <- params1$model.code
  axes <- c("X", "Y")
  
  for (ax in axes) {
    # Extract loadings
    loads1 <- model1$loadings[[ax]] %>%
      as.data.frame() %>% tibble::rownames_to_column(., "var")
    loads2 <- model2$loadings[[ax]] %>%
      as.data.frame() %>% tibble::rownames_to_column(., "var")
    
    plts <- list()
    # Iterate over components
    for (comp in 1:ncomp) {
      comp <- paste0("comp", comp)
      
      # Merge datasets for ease of use
      dat <- dplyr::inner_join(loads1 %>% dplyr::select(., c(var, comp)), 
                               loads2 %>% dplyr::select(., c(var, comp)), 
                               by = c("var")) %>%
        dplyr::mutate(size = sqrt(.[[2]]^2 + .[[3]]^2)) %>%
        dplyr::mutate(col = abs(.[[2]] - .[[3]]))
      
      title <- paste0("Loadings for ", ax, ".")
      subtitle <- paste0("Model: ", model.code)
      plt <- ggplot2::ggplot(data = dat, mapping = aes_string(x = paste0(comp, ".x"), 
                                                              y = paste0(comp, ".y"), 
                                                              label = "var", 
                                                              color = "col")) +
        ggplot2::geom_point(data = dat, mapping = aes(size = abs(5 * size))) +
        ggplot2::geom_abline(slope = 1, intercept = 0) +
        ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_vline(xintercept = 0) +
        ggrepel::geom_text_repel(data = dat) +
        ggplot2::xlab(paste0(comp, " - T1")) + ggplot2::ylab(paste0(comp, " - T2")) +
        ggplot2::labs(title = title, subtitle = subtitle)
      
      plts <- append(plts, list(plt))
    } # End loop components
    
    ret <- gridExtra::grid.arrange(grobs = plts, ncol = 2)
    ggplot2::ggsave(plot = ret, 
                    filename = paste0("loadings_", ax, "_", model.code, ".png"), 
                    device = "png", 
                    path = paste0("./results/compare_time/"), 
                    dpi = 270, height = 20, width = 20)
  } # End loop axes
  
}
