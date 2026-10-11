#' save_geoglist
#'
#' Save a `geoglist` in R's .RData format. As the raster or vector components of
#' a `geoglist` cannot be saved directly, this function handles restructuring of
#' that data to a format retrievable by the complementary function
#' `read_geoglist()`. These functions may be convenient to for storing and
#' retrieving large `geoglist` objects which took significant time to create.
#'
#' @param geog `geoglist`. A `geoglist` produced by `rast_to_geoglist()`,
#' optionally with island links.
#' @param path `character`. A character string specifying the directory and
#' file name to be used to store or retrieve the `geoglist`. This should end in
#' ".RData"
#' @return No return value.
#' @export

save_geoglist <- function(geog, path) {

  #geog <- foo
  #path <- "data/test.RData"

  if(!exists("geog")) {
    stop("Please supply geog as a geoglist")
  }
  if(!inherits(geog, "geoglist")) {
    stop("Please supply geog as a geoglist")
  }
  if(!exists("path")) {
    stop("Please supply path as a single character string")
  }
  if(!is.atomic(path)) {
    stop("Please supply path as a single character string")
  }
  if(!is.character(path) | length(path) != 1) {
    stop("Please supply path as a single character string")
  }
  if(!grepl("\\.RData$", path)) {
    stop("Path should end in '.RData'")
  }

  if(inherits(geog$layers, "SpatVectorCollection")) {
    geog$layers <- do.call(cbind, lapply(geog$layers, `[[`, 1))
    colnames(geog)
  } else {
    geog$layers <- as.data.frame(geog$layers)
  }

  if(!is.null(geog$links)) {
    geog$links <- list(as.data.frame(geog$links), geom(geog$links))
  }

  # save data
  saveRDS(geog, path)
}
