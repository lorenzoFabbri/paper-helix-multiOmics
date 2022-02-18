# Ugly script to check that when creating the networks, the correct
# partial correlation values are used (e.g., from a matrix to a df and/or 
# tidygraph object)

load("../data/intermediate_res_ggm/ggms.RData")
load("../data/intermediate_res_ggm/processed_ggms.RData")

time.point <- 2
pcor <- ggms[[time.point]]$model %>%
  as.matrix() %>%
  `class<-`("numeric") %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "mol") %>% tibble::as_tibble()
net <- processed.ggms[[time.point]]$net

node.attributes <- net %>% tidygraph::activate(., what = "nodes") %>%
  tidygraph::as_tibble() %>% dplyr::mutate(idx = seq.int(nrow(.)))
edge.attributes <- net %>% tidygraph::activate(., what = "edges") %>%
  tidygraph::as_tibble()
df <- dplyr::inner_join(node.attributes, edge.attributes, 
                        by = c("idx" = "from"))
df <- dplyr::inner_join(node.attributes, df, 
                        by = c("idx" = "to"))

if (dim(df)[1] != dim(edge.attributes)[1]) { stop(call. = TRUE) }

how.many <- 500
samples <- sample.int(dim(df)[1], how.many)
for (idx in samples) {
  node.a <- df[idx, ][["name.x"]]
  node.b <- df[idx, ][["name.y"]]
  val <- df[idx, ][["pcor"]] %>% as.numeric()
  
  to.check <- pcor %>%
    dplyr::filter(mol == node.a) %>%
    dplyr::select(node.b) %>% as.numeric()
  
  if (!(abs(val - to.check) < 1e-9)) {
    stop(call. = TRUE)
  }
}
