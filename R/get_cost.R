#' get_cost
#'
#' Get the cost of movement between cells in a `tardis` graph. Costs can be
#' calculated between pairs of cell, from one cell to many, or pairwise between
#' a set of cells.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param weights `character`. The name of the weighting scheme in `tardis` to
#' use for distance calculation. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of another weighting scheme added to
#' `tardis` using `weight_tardis()`.
#' @param origin `character`. A vector of `tardis` cell IDs.
#' @param dest `character`. A vector of `tardis` cell IDs.
#' @return For a pair of cells, a single value is returned. For one cell to
#' many, a vector of cost values is returned. For a set of cells, a square
#' matrix of costs is returned
#' @import cppRouting
#' @export
#'
#' @details
#' Getting the costs for a set of cells effectively creates a distance matrix.
#' This can very quickly create large objects. For example, a `tardis` graph
#' with ~20,000 accessible cells will give a 20k x 20k matrix approaching a
#' gigabyte in size. You may wish to optimise the resolution and masking of
#' your landscapes first if costing a large set of cells.
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
#' pairs <- get_cost(htd, origin = "19", dest = "1806")
#' one_to_many <- get_cost(htd, origin = "19", dest = c("19", "23", "26", "31"))
#' many_to_many <- get_cost(htd, origin = c("19", "23", "26", "31"), dest = c("19", "23", "26", "31"))
#'}

get_cost <- function(tardis, weights = "gdist", origin, dest) {

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
  if(!is.atomic(origin)) {
    stop("origin should be a character vector of cell IDs")
  }
  if(!is.character(origin)) {
    stop("origin should be a character vector of cell IDs")
  }
  if(!all(origin %in% c(tardis$edges[,1:2]))) {
    stop("Not all elements of origin are valid cell IDs in tardis")
  }
  if(!is.atomic(dest)) {
    stop("dest should be a character vector of cell IDs")
  }
  if(!is.character(dest)) {
    stop("dest should be a character vector of cell IDs")
  }
  if(!all(dest %in% c(tardis$edges[,1:2]))) {
    stop("Not all elements of dest are valid cell IDs in tardis")
  }

  tardis <- instantiate_tardis(tardis = tardis, weights = weights)
  mat <- get_distance_matrix(tardis$tgraph, origin, dest)
  if(nrow(mat) == 1) {mat <- as.vector(mat[1,])}
  return(mat)
}
