# RSingleRust
This is an experimental package for devising R wrappers of SingleRust using extendr

See [https://github.com/SingleRust/SingleRust](https://github.com/SingleRust/SingleRust) for more information.

``` r
library(Seurat)
library(SingleCellExperiment)
library(RSingleRust)

# get sparse matrix from the SCE object
sce <- as.SingleCellExperiment(pbmc_small)
counts <- t(assay(sce, "counts"))

# get qc metrics from SingleRust
meta.data <- get_qc_metrics(counts,
                            "data",
                            paste0("cells", 1:nrow(counts)),
                            paste0("genes", 1:ncol(counts)))
meta.data <- as.data.frame(meta.data)
```
