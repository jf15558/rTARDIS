#' click_lcp
#'
#' Interactively calculate and display least cost paths by clicking their start
#' and end points on the landscape plotted by the function. Clicked points
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
#' @param n `integer`. the number of point pairs to run (only for click_lcp)
#' @param col `character`. The colour to use for plotting least cost paths.
#' @param ... Additional arguments passed to `plot.geoglist()`.
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
#' # click two start and end points on the map
#' #click_iso(tardis = htd, geog = hexes, time = 2, n = 2)
#' }

click_lcp <- function(tardis, weights = "gdist", geog, time = NULL, n = 1, col = "gold", ...) {

   #tardis = rtd
   #geog = rasts
   #time = 2
   #col = "gold"
   #n = 1

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
  dst <- cbind(click(n = n), rep(time, n))

  hpts <- point_check(tardis, rbind(org, dst))

  hlcp <- least_cost(tardis, weights = weights, origin = hpts[1:n,], dest = hpts[(n + 1):(n * 2),])

  plot(hpts, col = col, pch = 16, add = T)
  plot(hlcp, add = T, col = col, lwd = 2)
  #plot(st_wrap_dateline(hlcp$geometry, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")), add = T, col = col, lwd = 2)
}
