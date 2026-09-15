#' tardis_to_sparse
#'
#' Create either a sparse weighted adjacency or transition probability matrix
#' from a `tardis` graph.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param weights `character`. The name of the weighting scheme in `tardis` to
#' use for distance calculation. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of another weighting scheme added to
#' `tardis` using `weight_tardis()`.
#' @param mode `character`. One of `"adjacency"` or `"transition"`.
#' @return `matrix` A `Matrix::sparseMatrix` adjacency or transition probability
#' matrix.
#' @import cppRouting Matrix
#' @export
#'
#' @details
#' Thix function produces a matrix representation of the edgelist in `tardis`,
#' hence why it is sparse (see `get_cost()` for creating dense matrices). For
#' the adjacency matrix, the weights in the desired weighting scheme are taken
#' directly. For the probability matrix, however, the weights are normalised so
#' that they sum to one. As such, a huge weight for one edge versus a tiny
#' weight for another edge could end up with a similar proportional probability.
#'
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE), right = FALSE)
#'
#' hexes <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = TRUE, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' cs <- tardis_to_sparse(htd, mode = "transition")
#' aj <- tardis_to_sparse(htd, mode = "adjacency")
#'}

tardis_to_sparse <- function(tardis, weights = "gdist", mode = "adjacency") {

  #tardis <- h1
  #weights = "gdist"
  #mode = "hitting"

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }

  if(!is.atomic(weights) | length(weights) != 1) {
    stop("weights should only contain one element")
  }
  if(!is.character(weights)) {
    stop("weights should be a character string")
  }
  if(!weights %in% colnames(tardis$edges)) {
    stop("weights should be a column name in tardis$edges")
  }
  if(!is.atomic(mode) | length(mode) != 1) {
    stop("mode should only contain one element")
  }
  if(!is.character(mode) | !mode %in% c("adjacency", "transition")) {
    stop("mode should be one of 'adjacency' or 'transition'")
  }

  tardis <- instantiate_tardis(tardis = tardis, weights = weights)
  if(mode == "adjacency") {
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1,
                        x = tardis$tgraph$data$dist)
  } else {
    # matrix as conductance rather than resistance
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1, x = 1 / tardis$tgraph$data$dist)
    # normalise into probability matrix
    mat <- mat / rowSums(mat)
  }
  return(mat)
}
