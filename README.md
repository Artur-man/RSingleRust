# RSingleRust
This is an experimental package for devising R wrappers of SingleRust using extendr

See [https://github.com/SingleRust/SingleRust](https://github.com/SingleRust/SingleRust) for more information.

``` r
librart(RSingleRust)

# read dense matrix
read_matrix(matrix(runif(25),nrow = 5,ncol = 5),
            paste0("cells", 1:5),
            paste0("genes", 1:5))
```
