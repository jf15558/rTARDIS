#' rotation_path
#'
#' For a given starting cell and destination layer, get its position in all
#' intervening layers. Points can be shifted forwards or backwards in time
#' regardless of the time directionality in the `tardis` graph, but will only
#' produce meaningful results when locations are not homologous between layers
#' (i.e., where locations shift due to plate tectonic rotation or some
#' other process). If the cell location is lost in an intervening layer, the
#' path will terminate at its last available location, rather than in the
#' destination layer.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param points `numeric` or `SpatVector`. Either a numeric vector of cell IDs
#' a `SpatVector` of points with a cell ID variable named 'cell'.
#' @param layer `numeric`. The layer number to shift points to. If a single
#' number, then all points (regardless of their starting layer) will be shifted
#' to that layer. Alternatively a vector of numbers with as many elements as
#' points to permit different destination layers for each point.
#' @param time `numeric`. Instead of a layer number, numeric ages can be
#' supplied. These are resolved to layer numbers internally, then treated like
#' the `layer` argument.
#' @param as.lines `logical`. Should the point-to-point shifts in each layer be
#' returned as line objects? These may be convenient for plotting purposes.
#' `FALSE` by default.
#' @return A `SpatVector` of points or lines. Both record the input point from
#' which they are derived (`$feature`). For points, the corresponding cell ID 
#' (`$cell`) and graph layer (`$layer`) are recorded for each row-wise shift.
#' For lines, their start cell (`$from`), end (`$to`) cell, start layer (`$srt`)
#' and end layer (`$end`) are recorded instead.
#' @import sf terra h3jsr
#' @export
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE), right = FALSE)
#'
#' hexes <- rast_to_geoglist(gal, gal_m, times = c(seq(2.25, 0, -0.5), 0), as.hex = TRUE, hex = 6)
#' hexes <- link_islands(hexes)
#'
#' htd <- build_tardis(hexes)
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2))
#' hpts <- point_check(htd, org)
#'
#' hpts2 <- rotation_path(htd, hpts, time = 0)
#' }

rotation_path <- function(tardis, points, layer, time = NULL, as.lines = F) {

   #tardis = htd
   #points = hpts
   #time = 0
   #as.lines = T

  if(!exists("tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if(!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of build_tardis")
  }
  if(length(tardis$tdat) < 3) {
    stop("tardis only contains a single layer. No shifts available")
  }

  if(!exists("points")) {
    stop("Please supply points as a vector of TARDIS cell IDs")
  }
  if(inherits(points, "SpatVector")) {
    if(geomtype(points) != "points") {
      stop("points must be a SpatVector of points geometries only")
    }
    if(!"cell" %in% names(points)) {
      stop("If supplying points as a SpatVector, this must contain a column of numeric cell IDs named 'cell'")
    }
    if(!is.numeric(points$cell)) {
      stop("points$cell must be numeric")
    }
    if(any(is.na(points$cell))) {
      stop("points$cell cannot contain missing values")
    }
    pts <- points$cell
  } else {
    pts <- points
  }

  if(!is.atomic(pts) | !is.numeric(pts)) {
    stop("Please supply points as a vector of TARDIS cell IDs")
  }
  if(any(pts %% tardis$gdat[5] > tardis$gdat[5])) {
    stop("Some points exceed the maximum permitted cell layer number. Check if these points correspond to the supplied tardis object")
  }
  if(any(pts %/% tardis$gdat[5] > (length(tardis$tdat) - 1))) {
    stop("Some points exceed the number of layers in tardis. Check if these points correspond to the supplied tardis object")
  }
  tests <- as.character(pts) %in% tardis$tgraph$dict$ref
  if(!all(tests)) {
    stop(paste0(length(tests) - sum(tests), " points do not correspond to tardis cell IDs. Use point_check() to validate starting cells first"))
  }

  if(!exists("layer") & is.null(time)) {
    stop("One of layer or time must be supplied")
  }
  if(!is.null(time)) {
    if(!is.atomic(time)) {
      stop("time should be a singe positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(!is.numeric(time)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(any(time < 0)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(length(time) != 1 & length(time) != length(pts)) {
      stop("time should be a single positive number, or a vector with n target ages corresponding to n entries in points")
    }
    if(any(time < min(tardis$tdat) | time > max(tardis$tdat))) {
      stop("One or more elements in time does not fall within the temporal range of tardis")
    }
    layer <- sum(time < tardis$tdat)
  }

  if(!is.atomic(layer) | ! is.numeric(layer)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(length(layer) != 1 & length(layer) != length(pts)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(any(layer %% 1 != 0 | layer < 0)) {
    stop("layer should be a single positive integer, or a vector with n target layers corresponding to n entries in points")
  }
  if(any(layer > (length(tardis$tdat) - 1))) {
    stop("One or more elements in layer exceeds the number of layers in tardis")
  }
  if(length(layer) != length(pts)) {
    rep(layer, length(pts))
  }

  pts <- cbind(1:length(pts), pts, (pts %/% tardis$gdat[5]) + 1, layer)
  rot <- tardis$edges[which(tardis$edges[,3] == 2),1:2]
  if(!is.na(tardis$gdat[7])) {
    grd <- get_grid(tardis$gdat[1:4], tardis$gdat[7])
  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
  }

  pth <- list()
  for(i in 1:nrow(pts)) {

    pth1 <- pt <- pts[i,2]
    bns <- ((pts[i] %/% tardis$gdat[5]) + 1):pts[i,4]
    for(k in bns) {
      # forwards in time
      if(bns[1] < bns[2]) {
        rt <- rot[which(((rot[,1] %/% tardis$gdat[5]) + 1) == k),1:2]
        # backwards in time
      } else {
        rt <- rot[which(((rot[,2] %/% tardis$gdat[5]) + 1) == k),2:1]
      }
      pt <- rt[match(pt, rt[,1]),2]
      
      # if not NA, add to rotation path
      if(!is.na(pt)) {
        pth1 <- c(pth1, pt)
      }
    }
    if(length(pth1) == 1) {
      pth1 <- c(pth1, pth1)
    }
    if(!is.na(tardis$gdat[7])) {
      crd <- st_coordinates(cell_to_point(grd[pth1 %% tardis$gdat[5]]))[,1:2, drop = F]
    } else {
      crd <- xyFromCell(samprast, pth1)
    }
    if(as.lines) {
      geoms <- lapply(1:(nrow(crd) - 1), function(x) {
        vect(cbind(crd[x,], crd[x + 1,]), type = "lines")
      })
      geoms <- do.call(rbind, geoms)
      geoms$feature <- i
      geoms$from <- pth1[-length(geoms)]
      geoms$to <- pth1[-1]
      geoms$srt <- bns[1:length(geoms)]
      geoms$end <- bns[(1:length(geoms)) + 1]
      
    } else {
      geoms <- vect(crd)
      geoms$feature <- i
      geoms$cell <- pth1
      geoms$layer <- bns[1:length(geoms)]
    }
    pth[[i]] <- geoms
  }
  pth <- do.call(rbind, pth)
  return(pth)
}
