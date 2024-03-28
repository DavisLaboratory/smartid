## code to prepare `sim_sce_test` dataset goes here

### set seed for reproducibility
set.seed(123)
### set simulation parameters
sim_params <- splatter::newSplatParams(
  nGenes = 100,
  batchCells = 400,
  group.prob = rep(0.25, 4),
  de.prob = 0.1,
  # de.downProb = 0,
  de.facLoc = 0.5,
  de.facScale = 0.4
)

### simulate data
sim_sce_test <- splatter::splatSimulate(sim_params, method = "groups")
## get up markers based on fold change
fc <- 1
cols <- paste0("DEFacGroup", seq_along(unique(sim_sce_test$Group)))
defac <- as.data.frame(SummarizedExperiment::rowData(sim_sce_test)[, cols])
up <- lapply(cols, \(id)
             dplyr::filter(defac, if_all(-!!sym(id), \(x) !!sym(id)/x > fc)) |>
               rownames())
SummarizedExperiment::metadata(sim_sce_test)$up_markers <- setNames(up, cols)

### save data
usethis::use_data(sim_sce_test, overwrite = TRUE)
