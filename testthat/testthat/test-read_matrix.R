test_that("add assay", {

  # read dense matrix
  read_matrix(matrix(runif(25),nrow = 5,ncol = 5),
              paste0("cells", 1:5),
              paste0("genes", 1:5))
})
