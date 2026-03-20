# nolint start

#' @docType package
#' @usage NULL
#' @useDynLib RSingleRust, .registration = TRUE
NULL

#' Return string `"Hello world!"` to R.
#' @export
hello_world <- function() .Call(wrap__hello_world)

#' Return string `"Hello world!"` to R.
#' @export
read_h5ad_memory <- function(file_path) .Call(wrap__read_h5ad_memory, file_path)

# nolint end
