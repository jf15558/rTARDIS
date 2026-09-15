#' cretaceous
#'
#' Load a `terra` R package `SpatRaster` object containing two palaeogeographic
#' digital elevation models of the Earth during the Cretaceous. Models are in
#' geographic coordinates with elevations in metres relative to sea level.
#' @return `SpatRaster`.
#' @export
#'
#' @details Each DEM records the estimated topography and bathymetry of the
#' Earth at 115 and 110 million years in the past. As such, each time slice is
#' taken to represent the 2.5 million years either side of these points (i.e.,
#' the 110 Ma reconstruction spans 112.5-107.5 Ma). Both reconstructions come
#' from the PALEOMAP plate rotation model series by Chris Scotese and
#' colleagues, and were obtained using the `fetch()` function of the
#' `chronosphere` R package.
#'
#' @source Kocsis, A. Scotese, C., Raja, N. (2013) chronosphere: Evolving
#' Earth System Variables. Zenodo. https://doi.org/10.5281/zenodo.3525481.
#'
#' Scotese, C., Verard, C., Burgener, L., Elling, R., Kocsis, A. (2024).
#' Phanerozoic-scope supplementary material to 'The Cretaceous World: Plate
#' Tectonics, Paleogeography, and Paleoclimate' from the PALEOMAP project
#' (Version v24221). Zenodo. https://doi.org/10.5281/zenodo.10659112.
#'
#' Scotese, C., Verard, C., Burgener, L., Elling, R., Kocsis, A. (2025). The
#' Cretaceous World: Plate Tectonics, Paleogeography, and Paleoclimate.
#' Geological Society, London, Special Publications, 544, SP544b

cretaceous <- function() {
  rast(system.file("extdata", "cretaceous.tif", package = "rTARDIS"))
}
