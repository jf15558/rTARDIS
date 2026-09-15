#' slice_geoglist
#'
#' Subset a `geoglist` to a given set of layers, along with any links if
#' present. If temporal. A discontiguous set of layers can be selected, unlike
#' with `slice_tardis()`, but any temporal information will be discarded if
#' so.
#'
#' @param geog `geoglist`. The output of `rast_to_geoglist()`.
#' @param times `numeric`. A vector of numbers denoting the desired times of
#' layers to subset from `geog`. One of `times` or `layers` must be specified.
#' Times corresponding to duplicated layer numbers will be ignored.
#' @param layers `numeric`. A vector of numbers denoting the layers to subset
#' from `geog`. One of `times` or `layers` must be specified. Duplicated layer
#' numbers will be ignored.
#' @return A `geoglist` of the requested subset of layers.
#' @export
#'
#' @examples
#' \donttest{
#' library(terra)
#' library(rTARDIS)
#'
#' gal <- galapagos()
#' gal_m <- classify(gal, matrix(c(-Inf, 0, NA, 0, Inf, 1), ncol = 3, byrow = TRUE), right = FALSE)

#' rasts <- rast_to_geoglist(gal, gal_m, times = c(seq(2.25, 0, -0.5), 0))
#' rasts1 <- slice_geoglist(rasts, layers = 1)
#' rasts2 <- slice_geoglist(rasts, layers = c(1, 3))
#' }

slice_geoglist <- function(geog, times = NULL, layers = NULL) {

  if(!inherits(geog, "geoglist")) {
    stop("geog should be `geoglist` object")
  }

  # check subsetting conflict
  if(is.null(times) & is.null(layers)) {
    stop("One of times or layers must be not be NULL")
  }
  if(!is.null(times) & !is.null(layers)) {
    stop("One of times or layers must be left as NULL")
  }

  # check times
  if(!is.null(times)) {

    if(is.null(geog$tdat)) {
      stop("geog contains no temporal information. Use the layers argument to subset instead")
    }
    if(!is.numeric(times)) {
      stop("Please supply times as a vector of numbers denoting the desired layer times to extract from geog")
    }
    if(any(is.na(times))) {
      stop("times cannot contain NA values")
    }
    if(any(times < 0)) {
      stop("times must only contain values >= 0")
    }
    times <- times[order(times, decreasing = T)]
    if(times[1] > geog$tdat[1] | times[length(layers)] < geog$tdat[length(geog$tdat)]) {
      stop("times must fall within the temporal range of geog")
    }
    layers <- c(sum(times[1] <= geog$tdat), sum(times[2] < geog$tdat))
  }

  # check layers
  if(!is.null(layers)) {

    if(!is.numeric(layers)) {
      stop("Please supply layers as an integer vector of the desired geoglist layers")
    }
    if(any(is.na(layers))) {
      stop("layers cannot contain NA values")
    }
    if(any(layers < 0) | any(layers %% 1 != 0)) {
      stop("layers must only contain positive integers")
    }
    layers <- layers[order(layers, decreasing = F)]
    if(layers[1] > length(geog$tdat) - 1 | layers[length(layers)] > length(geog$tdat) - 1) {
      stop("The values in layers cannot exceed the number of layers in geog")
    }
  }

  if(!is.atomic(layers)) {
    stop("layers should be a numeric vector of layer(s) to extract")
  }
  if(!is.numeric(layers)) {
    stop("layers should be a numeric vector of layer(s) to extract")
  }
  if(any(layers %% 1 != 0)) {
    stop("All elements of layers should be positive integers")
  }
  if(any(layers < 1)) {
    stop("All elements of layers should be positive integers")
  }
  nly <- ifelse(inherits(geog$layers, "SpatRaster"), nlyr(geog$layers), length(geog$layers))
  if(any(layers > nly)) {
    stop("One or more elements in layers exceeds the number of available layers in geog")
  }
  geog$layers <- geog$layers[unique(layers)]
  if(length(layers) == 1) {
    geog$layers <- svc(geog$layers)
  }
  if(!is.null(geog$links)) {
    geog$links <- geog$links[which(geog$links$layer %in% unique(layers))]
  }
  if(layers[length(layers)] - layers[1] != length(layers) - 1){
    warning("Requested layer set is not temporally contiguous. Temporal information will be discarded")
    geog$tdat <- NULL
  }
  return(geog)
}
