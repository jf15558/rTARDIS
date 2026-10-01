#' tardis_to_sparse
#'
#' Create either a sparse weighted adjacency or transition probability matrix
#' from a layer of a `tardis` graph.
#'
#' @param tardis `tardis`. The output of `build_tardis()`. This can only contain
#' a single time slice, so `slice_tardis()` may need to be used first.
#' @param weights `character`. The name of the weighting scheme in `tardis` to
#' use for distance calculation. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of another weighting scheme added to
#' `tardis` using `weight_tardis()`.
#' @param mode `character`. One of `"adjacency"` or `"transition"`.
#' @param raw.weight `logical`. If creating a transition matrix, should island
#' link probabilities be calculated from their raw weights (default) or instead
#' adjusted to approximate a jump-diffusion process (see details).
#' @param jump.dir `numeric` or `"auto"`. The stepwise cardinality if
#' approximating jump-diffusion probabilities. By default `"auto"`, in which
#' case the cardinality is selected based on the `tardis` landscape type.
#' Otherwise a suitable integer. Only relevant when `mode = "transition` and
#' `raw.weight = F` (see details for guidance).
#' @return `matrix` A `Matrix::sparseMatrix` adjacency or transition probability
#' matrix.
#' @import cppRouting Matrix
#' @export
#'
#' @details
#' This function produces a matrix representation of the edgelist in `tardis`,
#' hence why it is sparse (see `get_cost()` for creating dense matrices). For
#' the adjacency matrix, the weights in the desired weighting scheme are taken
#' directly.
#'
#' For the probability matrix, by default, the reciprocals of weights are
#' normalised so that they sum to one. As such, a huge weight for one edge
#' versus a tiny weight for another, unrelated edge could end up with similar
#' proportional probabilities. This behaviour is important for island links,
#' which will naturally be longer than links between adjacent cells.
#'
#' Island link probabilities will be proportionally smaller than other
#' transition probabilities from their starting cell, but can often display
#' higher probabilities for that single jump due to stepwise multiplication of
#' probabilities along alternative routes through unmasked space. This can be
#' accounted for by setting `jump = T`. Instead of using the proportional weight
#' directly, that weight is instead expressed as a multiple of the mean weights
#' for non-jump links from the focal cell (e.g., a link is 8.3 times the mean
#' non-link weight). A naive transition probability is then raised to the power
#' of this multiplier (m). By default this transition probability is selected
#' based on the structure of `tardis`. For a `tardis` built from a rectangular
#' raster grid, each grid cell has 8 neighbours, so the naive cost of moving in
#' any direction is 1/8 (= 0.125), raised to the power of m equivalent steps
#' along the link. For a hexagonally resampled `tardis`, unmasked cells have 6
#' neighbours. In this way, movement along the link is approximated as though it
#' were a series of direct steps through continuous space, where
#' other movement directions are probable but unrealised.
#'
#' 8 and 6 are sensible defaults for rectangular and hexagonal grids
#' respectively, but other values may be justified. For example `jump.dir = 2`
#' would approximate a scenario where a particle may move forwards or backwards
#' along a link with equal probability, but not deviate from the line of the
#' link itself, unlike with higher cardinalities.
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3,
#'                   byrow = TRUE), right = FALSE)
#'
#' hexes <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = TRUE, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' tr <- tardis_to_sparse(htd, mode = "transition")
#' aj <- tardis_to_sparse(htd, mode = "adjacency")
#' }

tardis_to_sparse <- function(tardis, weights = "gdist", mode = "adjacency",
                             raw.weight = TRUE, jump.dir = "auto") {

  #tardis <- htd
  #weights = "gdist"
  #mode = "transition"
  #jump = F
  #jump.dir = NULL

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.null(tardis$tdat)) {
    if(length(tardis$tdat) != 2) {
      stop("Only single-layer tardis objects can be converted to sparse. Use slice_tardis() to subset first")
    }
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

  if(!is.logical(raw.weights) | length(raw.weights) != 1) {
    stop("raw.weights should be a single logical value")
  }

  if(jump.dir != "auto") {
    if(!jump.dir %% 1 != 0) {
      stop("If not 'auto', then hex should be a an integer denoting the jump direction cardinality")
    }
  } else {

    cat(paste0("Autoselecting jump.dir from grid type"))
    jump.dir <- 8
    if(!is.na(tardis$gdat[7])) {jump.dir <- 6}
  }

  tardis <- instantiate_tardis(tardis = tardis, weights = weights)

  if(mode == "adjacency") {
    mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                        j = tardis$tgraph$data$to + 1,
                        x = tardis$tgraph$data$dist)
  } else {


    if(raw.weight) {

      # matrix as conductance rather than resistance
      mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                          j = tardis$tgraph$data$to + 1,
                          x = 1 / tardis$tgraph$data$dist)

      # normalise reciprocal raw weights into probability matrix
      mat <- mat / rowSums(mat)

    } else {

      # get edgelist, ordered by starting cell
      ed <- tardis$edges
      ed <- cbind(ed, 1:nrow(ed))
      ed <- ed[order(ed[,1]),]

      # get mean adjacent weight for each starting cell
      mean_wt <- ed[,weights]
      mean_wt[which(ed[,3] == 1)] <- NA
      mean_wt <- tapply(mean_wt, ed[,1], mean, na.rm = T)

      # if the edge type sum equals number of edges, it is a single cell without
      # adjacent weights to use, so take the global mean adjacent weight instead
      nedge <- table(ed[,1])
      singleton <- which(tapply(ed[,3], ed[,1], sum) == nedge)
      if(length(singleton) != 0) {
        mean_wt[singleton] <- mean(ed[which(ed[,3] == 0),weights])
      }

      # step multiplier for weight probability (adjacent weights unmodified)
      stepmult <- ed[,weights] / rep(mean_wt, nedge)
      stepmult[which(round(stepmult) == 1)] <- 1

      # cardinal probability compounded by the step multiplier gives jump prob
      jump_probs <- (1 / jump.dir) ^ stepmult

      # get remaining probability to distribute among adjacent cells:
      # (1 - sum of jump -probs, adjacent links having been zeroed)
      jump_probs[jump_probs == 1 / jump.dir] <- 0
      remainder <- tapply(jump_probs, ed[,1], function(x) {1 - sum(x)})

      # get proportional weights of adjacent cells (jumps having been zeroed)
      prop_wt <- ed[,weights]
      prop_wt[ed[,3] == 1] <- 0
      prop_wt <- unlist(tapply(prop_wt, ed[,1], function(x) {x / sum(x)}))

      # distribute remaining probability according to proportional weights
      # (jumps have zero probability because zero proportion of remainder)
      final <- rep(remainder, nedge) * prop_wt

      # restore jump probabilities to the values calculated earlier
      final[final == 0] <- jump_probs[final == 0]

      # validate all sets sum to 1
      if(!all(tapply(final, ed[,1], sum) == 1)) {
        stop("Error during jump probability calculation - there is a bug to squash")
      }

      # create sparse matrix
      final <- final[order(ed[,ncol(ed)])]
      mat <- sparseMatrix(i = tardis$tgraph$data$from + 1,
                          j = tardis$tgraph$data$to + 1,
                          x = final)
    }
  }
  return(mat)
}
