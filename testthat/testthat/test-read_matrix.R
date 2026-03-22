test_that("read dense matrix", {

  # read dense matrix
  read_matrix(matrix(runif(25),nrow = 5,ncol = 5),
              paste0("cells", 1:5),
              paste0("genes", 1:5))
})

test_that("read ondisk csr matrix (hdf5)", {

  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("anndataR")

  # write sparse matrix to h5
  mat <- matrix(sample(c(0,1), 25, replace = TRUE), nrow = 5, ncol = 5)
  mat <- as(mat, "dgRMatrix")
  tmp <- tempfile(fileext = ".h5")
  file <- rhdf5::H5Fcreate(tmp)
  anndataR:::write_h5ad_sparse_array(mat, file = file, name = "X", compression = "none")
  rhdf5::h5closeAll()

  # read sparse
  read_matrix_hdf5(tmp, "X/")
})

test_that("read ondisk csr matrix (anndata)", {

  skip_if_not_installed("zellkonverter")

  # zellkonverter
  h5ad_file <- system.file("extdata", "example_anndata.h5ad",
                           package="zellkonverter")

  # read sparse matrix
  read_matrix_hdf5(h5ad_file, "/obsp/connectivities")
})

test_that("read ondisk csr matrix (large hdf5)", {

  skip_if_not_installed("HDF5Array")
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("anndataR")

  # write sparse matrix to h5
  mat <- rsparsematrix(1000000, 1000, density=0.1, repr = "R")
  tmp <- tempfile(fileext = ".h5")
  file <- rhdf5::H5Fcreate(tmp)
  anndataR:::write_h5ad_sparse_array(mat, file = file, name = "X", compression = "none")
  rhdf5::h5closeAll()

  # read sparse
  read_matrix_hdf5(tmp, "X/")
})
