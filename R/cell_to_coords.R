#' cell_to_coords
#'
#' Convert a cell integer ID to its lon-lat coordinates, spatial cell ID and
#' layer, based on based on the general gridding parameters used by `rTARDIS`.
#' This will work for almost any integer (see details) and does not relate to
#' any specific landscape structure. See `point_check()` instead if you are
#' trying to resolve numeric coordinates to cells within a specific landscape.
#'
#' @param cellid `vector`. A vector of cell integer IDs, either numeric or
#' character.
#' @param gdat Either a vector in gdat format (see examples), or an object
#' containing this information (i.e., `tardis` or `geoglist`).
#' @param modulo.warning `logical`. Should the function warn the user if zero
#' modulo cells are present?
#' @return A `data.frame` recording the lon-lat coordinates, spatial cell
#' IDs and layer numbers of the input cell IDs. Integers which modulo to zero
#' will return `NA` for their lon-lat coordinates (see details).
#' @importFrom terra rast
#' @importFrom terra xyFromCell
#' @importFrom h3jsr cell_to_point
#' @export
#'
#' @details
#' Spatial cell position is retrieved as the modulo of layer cell count:
#' `cellid %% gdat[5]`. Temporal position is retrieved by integer division of
#' the layer cell count: `(cellid %/% gdat[5]) + 1`.
#'
#' In very rare cases (one in tens of thousands of cell), the modulo operation
#' results in a spatial cell index of zero and will not work. This behavior is
#' accounted for in all `rTARDIS` functions by masking such cells from the
#' landscape. However, lon-lat coordinates for these rare cells cannot be
#' derived, hence why they are returned as `NA` in the `data.frame` output.
#'
#' @examples
#' library(terra)
#' library(rTARDIS)
#'
#' # using an rTARDIS object
#' gal <- galapagos()
#' rasts <- rast_to_geoglist(gal[[1]])
#' cell_to_coords(1012, rasts)
#'
#' # using a manually-built gdat
#' gvec <- setNames(
#'              c(-92.0048120, -88.0048120, -2.0010895, 0.9989105, 3609, NA, 6),
#'              c("xmin", "xmax", "ymin", "ymax", "ncell", "ncol", "hex")
#'              )
#' cell_to_coords(c(1, 5, 233, 18), gvec)

cell_to_coords <- function(cellid, gdat, modulo.warning = TRUE) {

  cellid <- as.numeric(cellid)
  if(any(is.na(cellid))) {
    stop("Some cell IDs are not, or cannot be coerced to numeric")
  }
  if(!all(cellid %% 1 == 0)) {
    stop("Some cell IDs are not integers")
  }

  if(!inherits(gdat, "tardis") & !inherits(gdat, "geoglist")) {
    if(!is.vector(gdat)) {
      stop("gdat is not a tardis or geoglist, nor conforms to a $gdat vector")
    }
    if(length(gdat) != 7 | !all(names(gdat) == c("xmin", "xmax", "ymin", "ymax", "ncell", "ncol", "hex"))) {
      stop("gdat is not a tardis or geoglist, nor conforms to a $gdat vector")
    }
  } else {
    gdat <- gdat$gdat
  }

  pos <- cellid %% gdat[5]
  pos[pos == 0] <- NA
  lyr <- (cellid %/% gdat[5]) + 1

  if (!is.na(gdat[7])) {
    grid <- get_grid(gdat[1:4], gdat[7])
    crd <- st_coordinates(cell_to_point(grid[pos], gdat[7]))
    crd[is.na(pos),] <- NA

  } else {
    samprast <- rast(nrows = gdat[5] / gdat[6], ncols = gdat[6], ext = ext(gdat[1:4]))
    crd <- xyFromCell(samprast, pos)
    crd[is.nan(crd[,1])] <- NA
  }

  if(any(is.na(crd))) {
    if(modulo.warning) {
      warning("Some cells modulo to zero and will return NA coordinates")
    }
  }

  res <- cbind.data.frame(crd, pos, lyr)
  colnames(res) <- c("lon", "lat", "cell", "layer")
  return(res)
}
