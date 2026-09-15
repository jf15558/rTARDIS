#' rast_to_geoglist
#'
#' Convert a set of rasters to an S3 class `geoglist` object compatible with
#' downstream `rTARDIS` functions. Typically, these rasters will record
#' topography and/or bathymetry measured in metres, but could record other
#' geographic properties instead. If the rasters represent the same area through
#' time, then the first raster in the set must be the oldest.
#'
#' @param geog `SpatRaster`. A set of geographic rasters. These must be in
#' longitude-latitude projection. Missing values are not permitted and should
#' be replaced (e.g. with the average of the surrounding cells, or a dummy
#' value).
#' @param mask `SpatRaster` or `NULL`. If not `NULL`, this will be used to
#' designate non-accessible areas in `geog`. It must be fully contiguous with
#' `geog` (i.e., share the same resolution, extent and number of layers) and
#' contain only `1` (non-masked, accessible) or `NA` (masked, non-accessible)
#' values.
#' @param times `numeric` or `NULL`. If the layers do not relate to each other,
#' temporally, or only a single raster layer is present, then `times` can be
#' left as `NULL`. Otherwise a numeric vector with `nlayers(geog) + 1` positive
#' elements, with successive pairs expressing the temporal extent of each layer
#' as time in the past, such that: `t(n) > x >= t(n + 1)`. For temporal
#' ordering, the vector need not end in the present (i.e. `0`), but time must
#' flow from oldest to youngest.
#' @param as.hex `logical`. Should `geog` be resampled to the hexagonal grid
#' system defined by Uber's H3 library? Defaults to `FALSE`.
#' @param hex `"auto"` or `integer`. The desired H3 resolution to be used for
#' resampling rasters. Defaults to `"auto"`, resulting in the function selecting
#' the H3 resolution closest to the resolution of `geog`. Otherwise an integer
#' in the range 1 - 15.
#' @param method `character`. The function to be used for resampling the raster
#' grid. This must be compatible with `exactextractr::exact_extract()`.
#' Defaults to `"mean"`.
#' @param verbose `logical`. Should function progress be reported to the user?
#' @param ... Additional arguments passed internally to
#' `exactextractr::exact_extract()` for resampling of raster grids.
#' @return A `geoglist` with four list elements. `$gdat` records spatial
#' properties of the input rasters used throughout downstream TARDIS functions.
#' `$tdat` records the temporal extent of each layer, or is `NULL` if this
#' information was omitted in the function run. `$layers` is a set of geographic
#' layers, either as a standard `SpatRaster`, or a `SpatVectorCollection` of
#' hexagonal polygons if resampling was implemented. `$links` is a `NULL`
#' placeholder slot for storing the output of linking functions used later on.
#' @import terra exactextractr h3jsr
#' @export
#'
#' @details
#' Masking is a key feature of `rTARDIS`. Besides providing more realistic
#' representation of accessible geographic space, it can dramatically reduce the
#' number of cells and so memory required for landscape representation,
#' improving computational efficiency.
#'
#' Resampling to a hexagonal grid is implemented using
#' `exactextractr::exact_extract()` with weighting by cell area, as this is
#' currently much faster compared to the weighted resampling in the `terra` R
#' package. Weighting is not meaningful for all functions, e.g., `"max"`.
#'
#' Resampling to a hexagonal grid may be desirable when working with landscapes
#' approaching global extents in reduce the number of cells required to
#' represent polar latitudes, again helping to improve computational efficiency.
#' The trade off is that landscape features may be altered or lost depending on
#' the grid resolution used, although this a risk of any resampling procedure.
#' In addition, a large rectangular grid layer may be noticeably faster to plot
#' than a hexagonal grid of similar extent and resolution due to their handling
#' by `terra`, affecting performance in some downstream functions.
#'
#' @examples
#' \donttest{
#' # load libraries
#' library(terra)
#' library(rTARDIS)
#'
#' # load a dataset of the Galapagos archipelago through geological time
#' gal <- galapagos()
#'
#' # create a land-sea mask from the archipelago raster set
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE), right = FALSE)
#'
#' # create a geoglist from a single raster layer
#' rasts <- rast_to_geoglist(gal[[1]], gal_m[[1]])
#'
#' # create a multi-layer geoglist with hexagonal resampling
#' hexes <- rast_to_geoglist(gal, gal_m, times = c(seq(2.25, 0, -0.5), 0), as.hex = TRUE, hex = 6)
#' }

rast_to_geoglist <- function(geog, mask = NULL, times = NULL, as.hex = FALSE, hex = "auto", method = "mean", verbose = TRUE, ...) {

  #gal <- galapagos()
  #gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)

  #geog  = cret
  #mask = cret_l
  #as.hex = T
  #times = NULL
  #hex = 2
  #method = "mean"
  #verbose = T

  # check geography
  if(!exists("geog")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!inherits(geog, "SpatRaster")) {
    stop("Supply geog as a SpatRaster")
  }
  if(!is.lonlat(geog)) {
    stop("geog should be in geographic (long-lat) projection")
  }
  crs(geog) <- "EPSG:4326"
  if(any(is.na(geog[]))) {
    stop("NA values present in geog")
  }

  # check temporal data if supplied
  if(!is.null(times)) {
    if (!is.numeric(times) | length(times) != nlyr(geog) + 1) {
      stop("times must be a vector of time bin boundaries with nlayers(geog) + 1 elements")
    }
    if (any(diff(times) > 0)) {
      stop("All elements of times should be positive (i.e. before present) and in descending age order")
    }
  }

  if(!is.null(mask)) {

    if (!inherits(mask, "SpatRaster")) {
      stop("If not NULL, supply mask as a SpatRaster")
    }
    if (any(!unique(mask[]) %in% c(1, NA))) {
      stop("Mask layers can only contain 1 or NA values")
    }
    if (!is.lonlat(mask)) {
      stop("mask should be in geographic (lon-lat) projection")
    }
    crs(mask) <- "EPSG:4326"
    if (!all(dim(geog) == dim(mask)) | ext(geog) != ext(mask)) {
      stop("geog and mask must all have the same extent, resolution, and number of layers")
    }
    geog <- mask(geog, mask)
  }
  if(!is.atomic(as.hex) | length(as.hex) != 1) {
    stop("as.hex should be a logical indicating if to resample the raster grid to a hexagonal grid")
  }

  if(as.hex) {

    res <- max(cellSize(geog)[])
    res <- which.min(abs(res - h3jsr::h3_info_table$avg_area_sqm)) - 1
    if(!is.atomic(hex) | length(hex) != 1) {
      stop("hex should either be set to 'auto', or an integer specifying the desired H3 resolution for resampling [1-15]")
    }
    if(hex != "auto") {
      if(!hex %in% 0:15) {
        stop("If not 'auto', then hex should be a number specifying the desired H3 resolution for resampling [1-15]")
      }
      if(hex > res) {
        warning("The chosen H3 resolution is much smaller than the raster grid resolution. Resampling could take a while")
      }
    } else {
      cat(paste0("Autoselecting H3 resolution ", res, " "))
      hex <- res
    }

    if(!is.null(mask)) {

      mask <- disagg(mask, 2)
      grid <- lapply(mask, function(x) {xyFromCell(x, which(x[] == 1))})
    } else {
      mask <- disagg(geog[[1]])
      grid <- lapply(1:length(geog), function(x) {xyFromCell(mask, 1:ncell(mask))})
    }

    clist <- get_grid(as.vector(ext(geog)), hex)
    hex_list <- list()
    for(i in 1:nlyr(geog)) {

      if (verbose) {
        cat(paste0("Resampling layer [", i, "/", nlyr(geog), "]\r"))
      }

      # crop to prevent potential failures in retrieval of very high latitude cells
      cls <- unique(na.omit(unlist(suppressMessages(point_to_cell(grid[[i]], hex)))))

      # very rarely, some cells are recovered which do not lie in bounds - these are dropped
      cls <- intersect(cls, clist)
      clsp <- cell_to_polygon(cls)
      clsp <- st_make_valid(clsp)
      clsp <- st_wrap_dateline(clsp, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
      id <- match(cls, clist)
      vrs <- exact_extract(geog[[i]], clsp, fun = method, weights = "area")
      cat(paste0("\r"))
      dat <- data.frame(vrs)
      colnames(dat)[1] <- names(geog[[i]])
      dat$id <- id
      st_geometry(dat) <- clsp
      dat <- dat[order(dat$id),]

      ## CODE FOR MAKING SVC EQUIVALENT TO SPATRASTERSTACK
      foo <- vect(dat[st_is_valid(dat),])
      baz <- vect(cbind(1:(length(clist) - length(foo)), 1, NA, NA, 0), type = "polygon")
      baz$layer <- rep(NA, length(baz))
      baz$id <- setdiff(1:length(clist), foo$id)
      baz <- rbind(foo, baz)
      baz$layer[is.nan(baz$layer)] <- NA
      baz <- baz[order(baz$id)]
      baz$id <- clist

      hex_list[[i]] <- baz[,1]
    }
    out <- list(gdat = c(as.vector(ext(geog)), ncell = length(clist), ncol = NA, hex = hex), tdat = times, layers = svc(hex_list))

  } else {
    out <- list(gdat = c(as.vector(ext(geog)), ncell = ncell(geog), ncol = ncol(geog), hex = NA), tdat = times, layers = geog)
  }
  class(out) <- "geoglist"
  return(out)
}
