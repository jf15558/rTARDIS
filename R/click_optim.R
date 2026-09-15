#' click_optim
#'
#' Interactively calculate and display the optimal travel path linking a set of
#' points by clicking on the landscape plotted by the function. Clicked points
#' falling in masked regions are automatically resolved to the nearest available
#' cell.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param weights `character`. The name of the weighting scheme in `tardis` to
#' use for distance calculation. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of another weighting scheme added to
#' `tardis` using `weight_tardis()`.
#' @param geog `geoglist`. The `geoglist` used to build `tardis`.
#' @param time `integer`. The `tardis` time slice to plot and interact with.
#' Defaults to `NULL`, in which case the first slice is used.
#' @param n `integer`. The number of points to link with an optimal path. This
#' should be greater than 2.
#' @param loop `logical`. Should the optimal route additionally be closed from
#' end point to start point to return a polygon? Defaults to `FALSE`.
#' @param col `character`. The colour to use for plotting the optimal path.
#' @param ... Additional arguments passed to `plot.geoglist()`
#' @return No return value.
#' @import terra sf
#' @export
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE), right = FALSE)

#' hexes <- rast_to_geoglist(gal, gal_m, times = c(seq(2.25, 0, -0.5), 0), as.hex = TRUE, hex = 7)
#' hlink <- link_islands(hexes)
#' htd <- build_tardis(hexes)
#'
#' # click a point on the map
#' #click_iso(tardis = htd, geog = hexes, time = 2, cost = 1e5)
#' }

click_optim <- function(tardis, weights = "gdist", geog, time = NULL, n = 3, loop = FALSE, col = "gold", ...) {

  # tardis = rtd
  # geog = rasts
  # weights = "gdist"
  # time = 2
  # n = 3
  # col = "gold"
  # cost = 1e7

  if(n < 3) {
    stop("n should be greater than 2")
  }

  if(is.null(time)) {
    if(!is.null(tardis$tdat)) {
      time <- sum(tardis$tdat[1:2]) / 2
    }
    bin <- 1
  } else {
    if(time > tardis$tdat[1] | time < tardis$tdat[length(tardis$tdat)]) {
      stop("Time falls outside the range of tardis")
    }
    bin <- sum(time < tardis$tdat)
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  plot.geoglist(geog, bin, par.reset = F, ...)

  org <- cbind(click(n = n), rep(time, n))

  hpts <- point_check(tardis, org)

  hlcp <- optim_route(tardis, points = hpts, loop = loop)

  plot(hpts, col = col, pch = 16, add = T)
  plot(hlcp, add = T, col = col, lwd = 2)
}
