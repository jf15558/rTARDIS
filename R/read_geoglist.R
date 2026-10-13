#' read_geoglist
#'
#' Read a `geoglist` from a file produced by `save_geoglist()`. These functions
#' may be convenient to for storing and retrieving large `geoglist` objects
#' which took significant time to create.
#'
#' @param path `character`. A character string specifying the path to the
#' `geoglist` file. This should end in ".RData".
#' @return A `geoglist`, reconstructed from the data saved in `path`.
#' @import terra h3jsr
#' @export

read_geoglist <- function(path) {

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
  if(!file.exists(path)) {
    stop("Please supply path as a single character string")
  }

  # read stored geoglist
  geog <- readRDS(path)
  if(!inherits(geog, "geoglist")) {
    stop("path does not load an object of class geoglist")
  }

  # rebuild layers from data.frame
  if(is.na(geog$gdat[7])) {
    geog$layers <- rast(lapply(geog$layers, function(x) {
      rast(ncol = 360, nrow = 180, vals = x)
    }))

  } else {

    grd <- get_grid(geog$gdat[1:4], geog$gdat[7])
    not_na <- which(apply(geog$layers, 1, function(x) {any(!is.na(x))}))
    pols <- cell_to_polygon(grd[not_na])
    pols <- st_make_valid(pols)
    pols <- vect(st_wrap_dateline(pols, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")))

    pl <- lapply(1:ncol(geog$layers), function(x) {

      cat(paste0("Retrieving layers [", x, "/", ncol(geog$layers), "]\r"))
      if(x == ncol(geog$layers)) {
        cat("\n")
      }

      # build layer
      x1 <- geog$layers[[x]]
      present <- which(is.finite(x1))
      lyr <- pols[match(present, not_na)]
      lyr$val <- x1[present]
      lyr$id <- present

      # build blank geoms
      baz <- vect(cbind(1:(length(x1) - length(present)), 1, NA, NA, 0), type = "polygon")
      baz$val <- rep(NA, length(baz))
      baz$id <- setdiff(1:length(x1), lyr$id)

      # combine, order and return
      baz <- rbind(lyr, baz)
      baz <- baz[order(baz$id)]
      baz$val[is.nan(baz$layer)] <- NA
      baz$id <- NULL
      names(baz) <- names(geog$layers)[x]
      baz
    })
    geog$layers <- svc(pl)
  }

  # rebuild links if present
  if(!is.null(geog$links)) {
    lnk <- vect(geog$links[[2]], "lines")
    lnk$srt <- geog$links[[1]][,1]
    lnk$end <- geog$links[[1]][,2]
    lnk$layer <- geog$links[[1]][,3]
    lnk$distance <- geog$links[[1]][,4]
    geog$links <- lnk
  }

  # seems to help
  gc(verbose = F)

  # return geoglist
  return(geog)
}
