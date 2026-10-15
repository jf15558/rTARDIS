#' min_span
#'
#' Determine the least cost minimum spanning arborescence between sets of
#' points in a spatiotemporal landscape graph. Least cost paths are identified
#' between pairs of points using Djikstra's algorithm, then the minimum spanning
#' arborescence (see details) selected from these sets of pairwise distances.
#' Multiple sets of points may be analysed in one function call.
#'
#' @param tardis `tardis`. The output of `build_tardis()`.
#' @param weights `character`. The name of the weighting scheme in `tardis` to
#' use for distance calculation. By default these are true geographic distances
#' (`"gdist"`). Alternatively, the name of another weighting scheme added to
#' `tardis` using `weight_tardis()`.
#' @param points `SpatVector`. The output of `point_check()`, denoting the
#' points to be linked, and potentially with a grouping variable added
#' @param groups `character`. `NULL` by default, in which case all entries in
#' `points` will be linked in a single arborescence. Otherwise a column name in
#' `points` denoting subsets of points to link. This column can be `character`
#' or `numeric`, but must contain at least two points per group.
#' @param verbose `logical` Should function progress be reported to the user?
#' @return A `SpatVector` of time-discrete lines representing the branches of
#' the minimum spanning arborescences, recording their point group (`$group`),
#' the branch in the group-wise minimum spanning arboresence to which they
#' belong (`$feature`), the costs along each line (`$cost`) and their geographic
#' distances (`$distance`, identical if `weights = "gdist"`).
#' @import terra sf cppRouting rlemon h3jsr
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph as_edgelist
#' @importFrom igraph E
#' @importFrom igraph delete_edges
#' @importFrom rlemon MinCostArborescence
#' @export
#'
#' @details
#' If geographic distances are used (default), then there will be a single least
#' cost path between any two points and the minimum spanning arborescence will
#' be equivalent to a minimum spanning tree. However, as `rTARDIS` graphs are
#' directed, it is possible that the shortest path from A to B may not be the
#' same as the shortest path from B to A depending on the weighting scheme,
#' hence why minimum spanning arborescence is used as the more general case.
#' Note that the time directionality in `tardis` will also affect the minimum
#' spanning arborescence structure.
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
#' hmca <- min_span(htd, points = hpts)
#'
#' oldpar <- par(no.readonly = TRUE)
#' plot(hexes, par.reset = FALSE)
#' plot(hmca, add = TRUE, col = "gold", lwd = 2)
#' par(oldpar)
#' }

min_span <- function(tardis, weights = "gdist", points, groups = NULL, verbose = TRUE) {

  #tardis = ob
  #weights = "gdist"
  #points = pts
  #groups <- "tax"
  #verbose = TRUE

  if (!exists("tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if (!inherits(tardis, "tardis")) {
    stop("Supply tardis as the output of create_tardis")
  }
  if(!is.na(tardis$gdat[7])) {
    grid <- get_grid(tardis$gdat[1:4], tardis$gdat[7])

  } else {
    samprast <- rast(nrows = tardis$gdat[5] / tardis$gdat[6], ncols = tardis$gdat[6],
                     ext = ext(tardis$gdat[1:4]))
  }

  if(!is.atomic(weights) | length(weights) != 1) {
    stop("weights should only contain one element")
  }
  if(!is.character(weights)) {
    stop("weights should be a character string")
  }
  if(!weights %in% colnames(tardis$edges)) {
    stop("weights should be a column name in tardis$edges")
  }

  if(!class(points)[1] == c("SpatVector")) {
    stop("Supply origin as the output of point_check()")
  }

  if(is.null(groups)) {
    gps <- rep(1, length(points))
  } else {
    if(!is.atomic(groups) | length(groups) != 1) {
      stop("If not NULL, groups should be a single character string")
    }
    if(!is.character(groups)) {
      stop("If not NULL, groups should be a single character string")
    }
    if(!groups %in% names(points)) {
      stop("If not NULL, groups should be a column name in points")
    }
    if(any(table(points[[groups]][,1]) < 2)) {
      stop("If designating groups, each group must contain at least two points")
    }
    gps <- as.numeric(as.factor(points[[groups]][,1]))
  }

  if (!all(points$cell %in% tardis$tgraph$dict$ref)) {
    stop("One or more points in origin do not correspond to cells in tardis. Ensure that the correct tardis object is supplied for origin")
  }
  if (verbose) {
    cat("Initialising graph\n")
  }
  tardis <- instantiate_tardis(tardis = tardis, weights = weights)

  lst <- list()
  for(i in sort(unique(gps))) {

    if (verbose) {
      cat(paste0("Point group ", i, "/", max(gps), "\r"))
    }

    # distance matrix for points
    dists <- get_distance_matrix(tardis$tgraph, points$cell[gps == i], points$cell[gps == i])

    # set impossible paths to infinite and make get weighted subgraph for the point set
    dists[is.na(dists)] <- Inf
    gr <- graph_from_adjacency_matrix(dists, weighted = T)

    # drop inaccessible cells
    gr <- delete_edges(gr, which(E(gr)$weight == Inf))
    el <- as_edgelist(gr)
    el <- matrix(match(el, as.character(points$cell)), ncol = 2)
    ws <- E(gr)$weight

    # get the minimum spanning tree
    mca <- MinCostArborescence(el[,1], el[,2], ws, el[1,1], nrow(points))
    srt <- points$cell[mca$sources]
    end <- points$cell[mca$targets]

    # time points for cells
    srt_t <- srt %/% tardis$gdat[5] + 1
    end_t <- end %/% tardis$gdat[5] + 1

    # get pathways in minimum spanning tree
    paths <- lapply(suppressMessages(get_path_pair(tardis$tgraph, from = srt, to = end)),
                    function(x) {as.numeric(x)})

    costs <- lapply(1:length(paths), function(x) {
      if (verbose) {
        cat(paste0("Point group ", i, "/", max(gps), " [", ceiling((x / length(paths)) * 100), "%]   ", "\r"))
      }
      pth <- paths[[x]]
      wp <- c(get_distance_pair(tardis$tgraph, pth[-length(pth)],
                                pth[-1], algorithm = "Djikstra"), 0)
      if (weights != "tdist") {
        tp <- c(get_distance_pair(tardis$tgraph, pth[-length(pth)],
                                  pth[-1], algorithm = "Djikstra", aggregate_aux = T),
                0)
      }
      else {
        tp <- wp
      }
      list(wp, tp)
    })
    wpaths <- lapply(costs, function(x) {
      x[[1]]
    })
    tpaths <- lapply(costs, function(x) {
      x[[2]]
    })

    pords <- ifelse(srt_t - end_t <= 0, 1, 0)

    path_groups <- path_ids <- wvec <- tvec <- orgs <- dsts <- list()
    for (j in 1:length(paths)) {

      if(!is.na(tardis$gdat[7])) {
        path_xy <- st_coordinates(cell_to_point(grid[paths[[j]] %% tardis$gdat[5]]))

      } else {
        path_xy <- xyFromCell(samprast, paths[[j]] %% tardis$gdat[5])
      }

      pbin <- ceiling(paths[[j]] / tardis$gdat[5])
      pseq <- rep(1:length(rle(pbin)$length), rle(pbin)$length)
      wcst <- rev(rev(c(rbind(tapply(wpaths[[j]], INDEX = pseq,
                                     sum), 0)))[-1])
      dcst <- rev(rev(c(rbind(tapply(tpaths[[j]], INDEX = pseq,
                                     sum), 0)))[-1])
      pgrp <- list()
      within <- tapply(1:length(pseq), INDEX = pseq, function(y) {
        if (length(y) == 1) {
          st_linestring(path_xy[c(y, y), ])
        }
        else {
          st_linestring(path_xy[y, ])
        }
      })
      for (k in 1:length(seq(1, length(wcst), 2))) {
        pgrp[[seq(1, length(wcst), 2)[k]]] <- within[[k]]
      }
      pids <- paste0(j, "_", pbin[1], "-", pbin[1])
      if (length(unique(pbin)) > 1) {
        between <- lapply(which(diff(pseq) != 0) + 1, function(y) {
          st_linestring(as.matrix(path_xy[c(y - 1, y),
                                          , drop = FALSE]))
        })
        for (k in 1:length(seq(2, length(wcst), 2))) {
          pgrp[[seq(2, length(wcst), 2)[k]]] <- between[[k]]
        }
        pids <- pbin[c(1, which(diff(pbin) != 0) + 1)]
        pids <- rev(rev(c(t(cbind(paste0(j, "_", pids, "-",
                                         pids), c(paste0(j, "_", pids[-length(pids)],
                                                         "-", pids[-1]), NA)))))[-1])
      }
      pord <- if (pords[j] == 1) {
        1:length(pgrp)
      } else {
        rev(1:length(pgrp))
      }
      path_groups[[j]] <- pgrp[pord]
      path_ids[[j]] <- pids[pord]
      wvec[[j]] <- wcst[pord]
      tvec[[j]] <- dcst[pord]
    }
    path_ids <- unlist(path_ids)
    ob <- cbind.data.frame(group = unique(points[[groups]][,1])[i], feature = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                       "_"), function(y) {
                                                                         y[[1]]
                                                                       }))), srt_bin = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                                                                         "_|-"), function(y) {
                                                                                                                           y[[2]]
                                                                                                                         }))), end_bin = as.numeric(unlist(lapply(strsplit(path_ids,
                                                                                                                                                                           "_|-"), function(y) {
                                                                                                                                                                             y[[3]]
                                                                                                                                                                           }))), distance = unlist(tvec), cost = unlist(wvec))
    ob$layer <- (ob$srt_bin + ob$end_bin) / 2
    ob <- ob[,c("group", "feature", "layer", "cost", "distance")]
    st_geometry(ob) <- st_sfc(unlist(path_groups, recursive = F),
                              crs = "+proj=longlat")
    ob <- st_wrap_dateline(ob, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))

    if(verbose) {if(i == max(gps)) {cat("\n")}}
    lst[[i]] <- vect(ob)
  }

  return(do.call(rbind, lst))
}
