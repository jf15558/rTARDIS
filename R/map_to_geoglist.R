#' map_to_geoglist
#'
#' Convert a vector of values to a `geoglist`, using an existing `geoglist` as a
#' scaffold. This function is primarily intended to support mapping of outputs
#' from functions involving `tardis` objects and their derivatives, unlike
#' `rast_to_geoglist()` which creates their structure from scratch.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param values `numeric`. A named vector of numeric values. Missing values
#' like `NA` or `Inf` are permitted, although not necessarily encouraged.
#' Element names must be numeric IDs which correspond to non-masked cells from
#' a single layer in `geog`. Duplicate IDs are not permitted, i.e., only one
#' value to map per cell ID.
#' @param name `character`. The field name under which values will be stored.
#' @return A single layer `geoglist` with the mapped `values` present in
#' the `name` field.
#' @export
#'
#' @examples
#'
#' load galapagos dataset
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = T), right = F)
#'
#' # hex-gridded tardis
#' rasts <- rast_to_geoglist(gal[[1]], gal_m[[1]], as.hex = T, hex = 6)
#' rasts <- link_islands(rasts, klink = NULL)
#' rtd <- build_tardis(rasts)
#'
#' # small analysis
#' smc <- tardis_to_samc(tardis = rtd, absorption = 0.001)
#' dsp <- dispersal(smc, dest = smc@names[1], time = 1000)
#' names(dsp) <- foo@names
#'
#' foo <- map_to_geoglist(rasts, dsp)
#' plot(foo)

map_to_geoglist <- function(geog, values, name = "value") {

  geog = rasts
  values = dsp
  name = "value"

  if(!exists("geog")) {
    stop("Please supply geog as a geoglist")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Please supply geog as a geoglist")
  }

  if(!exists("values")) {
    stop("Please supply values as a named vector of numerics")
  }
  if(!is.atomic(values)) {
    stop("Please supply values as a named vector of numerics")
  }
  if(is.null(names(values)) | !is.numeric(values)) {
    stop("Please supply values as a named vector of numerics")
  }
  if(any(is.na(as.numeric(names(values))))) {
    stop("Some value names are not numbers and so cannot be resolved to geoglist cells")
  }

  if(!is.atomic(name) | length(name) != 1) {
    stop("name should be a single character string")
  }
  name <- as.character(name)

  cls <- as.numeric(names(values))
  if(any(duplicated(cls))) {
    stop("Duplicate cell IDs are not permitted")
  }
  if(any(cls %% geog$gdat[5] > geog$gdat[5])) {
    stop("Some value names exceed the highest permitted cell ID in geog")
  }

  nl <- length(geog$layers)
  if(is.na(geog$gdat[7])) {nl <- nlyr(geog$layers)}
  lyr <- unique((cls %/% geog$gdat[5]) + 1)
  if(length(lyr) != 1) {
    stop("Not all value names come from the same layer")
  }
  if(lyr > nl) {
    stop("Value names exceed the number of layers available in geog")
  }
  if(nl != 1) {
    geog <- slice_geoglist(geog, layers = lyr)
  }
  cls <- cls %% geog$gdat[5]

  if(any(is.na(geog$layers[[1]][[1]][cls,1]))) {
    stop("Some value names correspond to masked cells in geog")
  }

  elem <- geog$layers[[1]]
  elem[[1]] <- NA
  elem[[1]][cls,1] <- values
  names(elem) <- name

  if(inherits(geog$layers, "SpatRaster")) {
    geog$layers <- elem
  } else {
    geog$layers <- svc(elem)
  }

  return(geog)
}
