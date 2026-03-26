test_that("process single cell experiment", {

  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  # example
  sce <- as.SingleCellExperiment(pbmc_small)
  counts <- t(assay(sce, "counts"))

  # get qc metrics from SingleRust
  meta.data <- get_qc_metrics(counts,
                              "data",
                              paste0("cells", 1:nrow(counts)),
                              paste0("genes", 1:ncol(counts)))
  meta.data <- as.data.frame(meta.data)
})


