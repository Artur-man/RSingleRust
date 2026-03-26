test_that("process single cell experiment", {

  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")

  # example
  sce <- as.SingleCellExperiment(pbmc_small)
  # sparse_check_S4(sce@assays@data@listData$counts)
  # assay(sce, "counts") <- as.matrix(assay(sce, "counts"))

  # read
  counts <- t(assay(sce, "counts"))
  get_qc_metrics(counts,
                 "data",
                 paste0("cells", 1:nrow(counts)),
                 paste0("genes", 1:ncol(counts)))
})


