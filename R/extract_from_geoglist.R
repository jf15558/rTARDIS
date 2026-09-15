#' extract_from_geoglist
#'
#' Extract values from the layers in a `geoglist` using a terra R package
#' `SpatVector` containing one or several geometries. This will typically come
#' from another `rTARDIS` function (e.g, `least_cost()`, `isochrone()`, ect),
#' but could instead be user-designed.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param geom `SpatVector`. The geometry or geometries which will be used to
#' extract values from `geog`.
#' @param layer `numeric`. If not `NULL`, then an integer specifying the layer
#' in `geog` from which values are to be extracted. This argument is primarily
#' intended for use with user-designed `geom` objects which do not contain layer
#' assignments, unlike returns from other `rTARDIS` functions.
#' @return `SpatVector`. A `SpatVector` of points corresponding to the centroids
#' of all cells in `geog` intersected by each geometry in `geom` (denoted by
#' `$feature`), the layers of each intersected cell (`$layer`), and the cell
#' values (`$value`).
#' @import terra h3jsr
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
#' rasts <- rast_to_geoglist(gal, gal_m, times = c(seq(2.25, 0, -0.5), 0), as.hex = TRUE, hex = 6)
#' rasts <- link_islands(rasts)
#' rtd <- build_tardis(rasts)
#' org <- rbind(c(-89.78873, -1.420627, 2),
#'              c(-89.58525, -1.473917, 2))
#' dst <- rbind(c(-88.70836, -0.2627832, 2),
#'              c(-90.44276,  0.2943382, 2))
#'
#' rpts <- point_check(rtd, rbind(org, dst))
#' rlcp <- least_cost(rtd, origin = rpts[1,], dest = rpts[3,])
#' vals <- extract_from_geoglist(rasts, rlcp)
#' }

extract_from_geoglist <- function(geog, geom, layer = NULL) {

   #geog = rasts
   #geom = rlcp
   #layer = NULL

  if(!exists("geog")) {
    stop("Supply geog as a geoglist with rast_to_geoglist()")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Supply geog as a geoglist from rast_to_geoglist()")
  }
  if(!exists("geom")) {
    stop("Supply geom as an SpatVector object")
  }
  if(!is.null(layer)) {
    if(!is.atomic(layer) | length(layer) != 1) {
      stop("If not NULL, layer should be a single integer")
    }
    if(!is.numeric(layer)) {
      stop("If not NULL, layer should be a single integer")
    }
    if(!layer %% 1 != 0) {
      stop("If not NULL, layer should be a single integer")
    }
    geom$layer <- layer
  }
  if(is.null(geom$feature)) {
    geom$feature <- 1:nrow(geom$feature)
  }

  if(!inherits(geog$layers[[1]], "SpatRaster")) {
    grid <- get_grid(geog$gdat[1:4], geog$gdat[7])
  }

  vals <- lapply(1:max(geom$feature), function(x) {
    pth <- geom[which(geom$feature == x),]
    vals2 <- lapply(min(pth$layer):max(pth$layer), function(y) {
      prt <- pth[which(pth$layer == y),]
      lyr <- geog$layers[[y]]
      if(inherits(lyr, "SpatRaster")) {
        vl <- extract(lyr, vect(prt), cells = T)
        cls <- vect(xyFromCell(lyr, vl$cell))
        cls$feature <- rep(x, length(cls))
        cls$layer <- rep(y, length(cls))
        cls$value <- vl[,2]


      } else {

        vl <- which(relate(lyr, prt, "intersects")[,1])
        cls <- centroids(lyr[vl,1])
        cls$value <- cls[[names(cls)]][,1]
        cls$feature <- rep(x, length(cls))
        cls$layer <- rep(y, length(cls))
        cls <- cls[,c("feature", "layer", "value")]
      }
    })
    do.call(rbind, vals2)
  })
  return(do.call(rbind, vals))
}
