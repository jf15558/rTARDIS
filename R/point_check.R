#' point_check
#'
#' Resolve a set of coordinates to accessible cells in a `tardis` object.
#' Discrepant points are adjusted to their nearest accessible cell based on
#' great circle distance, with the option to cap the distance tolerated. Checks
#' are performed on unique coordinate sets internally for efficiency.
#'
#' @param tardis `tardis`. The output of `build_tardis()` or `weight_tardis()`.
#' @param points `matrix` or `data.frame`. A two or three column `matrix` or
#' `data.frame` of coordinates. Column ordering is assumed to be longitude
#' (decimal degrees), latitude (decimal degrees) and time (positive, time before
#' present). The time column is only required if `tardis` contains multiple
#' layers.
#' @param max.dist `numeric`. The maximum distance in metres permitted for
#' adjusting points to accessible cells. Points with adjustments above this
#' threshold will be discarded. `NULL` by default, meaning that no points are
#' discarded.
#' @param round `integer`. If not `NULL`, then the number of decimal places to
#' which longitudes and latitudes will be rounded. This can speed up function
#' progress when dealing with large numbers of input points.
#' @param verbose `logical` Should function progress be reported to the user?
#' @return A `SpatVector` of points, recording which input point they correspond
#' to in case points were discarded (`$feature`), the ID of the `tardis` cell
#' they fall into (`$cell`), the layer to which they belong (`$layer`) and their
#' adjusted distance (`$adj`).
#' @import geosphere terra sf h3jsr
#' @export
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
#' hexes <- rast_to_geoglist(gal[[1:2]], gal_m[[1:2]], times = c(2.25, 2, 1.75),
#'                           as.hex = TRUE, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#'
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2),
#'              c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' hpts <- point_check(htd, org)
#' }

point_check <- function(tardis, points, max.dist = NULL, round = NULL, verbose = TRUE) {

  #tardis <- ob
  #points <-  occs[,c("p_lng", "p_lat", "bin_midpoint")]
  #max.dist = NULL
  #round <- 2
  #verbose = TRUE

  if (!exists("tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if (!inherits(points, "data.frame") & !inherits(points, "matrix")) {
    stop("Supply origin as a data.frame or matrix with two or three columns")
  }
  if (ncol(points) < 2 | ncol(points) > 3) {
    stop("Supply points as a data.frame or matrix with two or three columns")
  }
  if(!is.null(max.dist)) {
    if(!is.atomic(max.dist) | length(max.dist) != 1) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
    if(!is.numeric(max.dist)) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
    if(max.dist <= 0) {
      stop("If not NULL, max.dist should be a single positive numeric")
    }
  }
  if(!is.null(round)) {
    if(!is.atomic(round) | length(round) != 1) {
      stop("If not NULL, round should be a single positive integer")
    }
    if(!is.numeric(round)) {
      stop("If not NULL, round should be a single positive integer")
    }
    if(round <= 0 | round %% 1 != 0) {
      stop("If not NULL, round should be a single positive round")
    }
  }

  if(length(tardis$tdat) < 3) {
    points <- cbind(points[,1:2], rep(0.5, nrow(points)))
    tardis$tdat <- c(1, 0)
  }
  if(any(is.na(points))) {
    stop("NA values are not permitted in points")
  }
  if(any(points[,1] > 180 | points[,1] < -180)) {
    stop("The first two columns in points must be in decimal degrees")
  }
  if(any(points[,2] > 90 | points[,2] < -90)) {
    stop("The first two columns in points must be in decimal degrees")
  }
  if(any(points[,3] < 0)) {
    stop("If supplied, the time column in points must only contain positive values")
  }
  if(length(tardis$tdat) > 2 & ncol(points) < 3) {
    stop("If tardis contains multiple time layers, then points should contain longitude, latitude and time columns")
  }
  if(is.null(tardis$tdat) & ncol(points) > 2) {
    stop("If tardis consists of one layer, then points should contain longitude and latitude columns only")
  }

  # trim to unique points, with rounding if selected
  if(!is.null(round)) {points[,1:2] <- round(points[,1:2], round)}
  id <- paste(match(points[,1], unique(points[,1])),
              match(points[,2], unique(points[,2])),
              match(points[,3], unique(points[,3])), sep = "_")
  id2 <- unique(id)
  points <- unique(points)

  # resolve coordinates to grid
  if(!is.na(tardis$gdat[7])) {
    grid <- get_grid(tardis$gdat[1:4], tardis$gdat[7])
    pcell <- suppressMessages(point_to_cell(points[, 1:2, drop = F], tardis$gdat[7]))
    pcell <- match(pcell, grid)

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
    pcell <- cellFromXY(samprast, points[, 1:2, drop = F])
  }
  pt <- sapply(points[, 3], function(y) {sum(y < tardis$tdat)})
  pt[which(points[, 3] > tardis$tdat[1] | points[, 3] < tardis$tdat[length(tardis$tdat)])] <- NA
  ptcell <- pcell + ((pt - 1) * tardis$gdat[5])
  pmod <- rep(0, length(ptcell))

  # get cells to resolve
  to_do <- which(!ptcell %in% tardis$edges[,1] & !is.na(ptcell))
  tms <- pt[to_do]

  # resolve per time slice
  counter <- 0
  for(i in unique(tms)) {

    # set to do
    pos <- to_do[which(tms == i)]
    if (verbose) {
      cat(paste0("Resolving discrepant points [", counter, "/", length(to_do), "]\r"))
      if (i == rev(unique(tms))[1]) {
        cat("\n")
      }
    }

    # get border cells
    cls <- tardis$edges[which(tardis$edges[,1] %/% tardis$gdat[5] == (i - 1) & tardis$edges[,3] == 0),1]
    cls <- table(cls)
    if(!is.null(tardis$gdat[7])) {
      cls <- names(cls[cls < 6])
    } else {
      cls <- names(cls[cls < 8])
    }

    # create point sets
    crd <- cell_to_coords(cls, tardis)
    crd <- vect(crd[,1:2], crs = "+proj=lonlat")
    crd2 <- vect(points[pos,1:2, drop = F], crs = "+proj=lonlat")

    # get nearest
    nr <- nearest(crd2, crd)
    dst <- nr$distance
    nr <- cls[nr$to_id]

    # add in new cells
    pmod[pos] <- dst
    ptcell[pos] <- nr
    counter <- counter + length(pos)
  }

  if(!is.null(max.dist)) {
    ptcell[pmod > max.dist] <- NA
  }

  if(all(is.na(ptcell))) {
    stop("No points could be matched to accessible cells in tardis")
  }
  if(any(as.logical(pmod))) {
    warning("Some points were translated to accessible cells in tardis")
  }
  if(any(is.na(pt))) {
    warning("Some points fell outside the temporal bounds of tardis and were discarded")
  }
  out <- cbind.data.frame(point = which(!is.na(ptcell)),
                          bin = pt[!is.na(ptcell)] %/% tardis$gdat[5] + 1,
                          mod = pmod[!is.na(ptcell)])

  if(!is.na(tardis$gdat[7])) {
    geom <- vect(cell_to_point(grid[pcell], tardis$gdat[7]))
    geom$cell <- match(geom$h3_address, grid)

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6], ext = ext(tardis$gdat[1:4]))
    geom <- vect(xyFromCell(samprast, pcell))
    geom$cell <- ptcell
  }
  geom$feature <- out$point
  geom$layer <- out$bin
  geom$adj <- out$mod

  # expand to original input and return
  geom <- geom[match(id, id2)]
  return(geom[,c(3,2,4,5)])
}
