# RSingleRust
This is an experimental package for devising R wrappers of SingleRust using extendr

See [https://github.com/SingleRust/SingleRust](https://github.com/SingleRust/SingleRust) for more information.

``` r
library(anndataR)
librart(RSingleRust)

file <- system.file("extdata", "example.h5ad", package = "anndataR")
read_h5ad_memory(file)
```
