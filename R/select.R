## This file is part of the rgl.cry package
##
## Functions

#' Select atoms or reciprocal lattice points.
#'
#' Select one or more atoms or reciprocal lattice points in the window.  The
#' labels and Miller indices of the selected atoms or lattice points will be
#' displayed.
#'
#' Selecting atoms or lattice points in the window will include all
#' z-coordinates.  If you do not want to include all z-coordinates, you will
#' need to modify the code.
#'
#' @param dev RGL device to apply.  Defaults to current device.
#' @param verbose logical: Should the report be suppressed?
#'
#' @return List of Miller indices or element labels.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'  select()
#'  select(dev = 123)
#' }
#' }
select <- function(dev = NULL, verbose = TRUE) {

  list(dev = dev, verbose = verbose)

  ## Select device
  if (missing(dev)) {
    tgt.dev <- rgl::cur3d()
  } else {
    tgt.dev <- dev
  }


  ## Get the current pair of demo and check if it exist.
  inst <- pkg$inst
  idx <- which(inst$dp.dev == tgt.dev | inst$cry.dev == tgt.dev)
  ## Since integer(0) is returned when the element is empty, it will be judged
  ## by length().

  if (length(idx) == 0) {
    stop("The device was lost.\n")
  }

  if (length(idx <- which(inst$dp.dev == tgt.dev)) != 0) {
    type <- "dp"
  } else if (length(idx <- which(inst$cry.dev == tgt.dev)) != 0) {
    type <- "cry"
  }


  ##
  if (verbose == TRUE) {
    message("To select points, drag the left mouse button.")
    message("To finish, press ESC.")
  }



  ## ------------------------------------------------------------
  ## Select reciprocal lattice points or atoms.

  cur.dev <- rgl::cur3d() # Save the current device.

  rgl::set3d(tgt.dev, silent = TRUE)

  if (type == "dp") {
    ## Limit objects to select.
    objIds <- rgl::ids3d(tag = TRUE, subscene = 0)
    objIds <- objIds[(objIds$tag == "rlpoints1"), ] # reciprocal lattice points.
    objIds <- objIds[(objIds$type == "spheres"), ] # just to be sure.
    ## If you limit the points based on visibility, ...

    selection <- rgl::selectpoints3d(objIds$id,
      value = FALSE,
      multiple = function(x) {
        ## x has id and index if `value' is
        ## FALSE, otherwise coordinate of points.
        ## cat(sprintf("id, idx: %0.f, %0.f\n",
        ##            x[,"id"], x[,"index"]))
        if (verbose == TRUE) {
          message(".", appendLF=FALSE)
        }
        ## spheres3d(x, radius = 0.05)
        TRUE
      }
    )
    if (verbose == TRUE) message()

    ## Reports
    hkl <- inst[[idx, "hkl"]]
    i <- selection[, "index"] # id, index (all ids are the same)
    i <- as.numeric(i)
    ret <- sapply(i, function(j) {
      sprintf("%s %s %s", hkl[j, "H"], hkl[j, "K"], hkl[j, "L"])
    })
  } else if (type == "cry") {
    ## Limit objects to select.
    objIds <- rgl::ids3d(tag = TRUE, subscene = 0)
    objIds <- objIds[(objIds$tag == "atom"), ] # atom
    objIds <- objIds[(objIds$type == "spheres"), ] # just to be sure.

    selection <- rgl::selectpoints3d(objIds$id,
      value = FALSE,
      multiple = function(x) {
        ## x has id and index if `value' is
        ## FALSE, otherwise coordinate of points.
        ## cat(sprintf("id, idx: %0.f, %0.f\n",
        ##            x[,"id"], x[,"index"]))
        if (verbose == TRUE) {
          message(".", appendLF=FALSE)
        }
        ## spheres3d(x, radius = 0.05)
        TRUE
      }
    )
    if (verbose == TRUE) message()

    rgl::set3d(cur.dev, silent = TRUE) # Restore the current device.

    ## The atoms are placed at pos[[i]], i ranges from 1 to the number of
    ## non-equivalent atoms.  The list containts the coordinates for each atom.
    ##
    ## The selection reports the id and index.  The id corresponds to the i
    ## above but not the same number.  The order in which they appear seems to
    ## be the same.  So it can be retrived as follows:

    ## Reports
    lCIF <- inst[[idx, "lCIF"]]
    ids <- unique(sort(objIds$id))
    i <- selection[, "id"] # id, index (id corresponds to atom)
    i <- as.numeric(i)
    ret <- sapply(i, function(j) {
      sprintf("%s", lCIF$COOR$VAL[which(ids == j), "label"])
    })

    ## print(objIds)
    ## print(selection)
    ## print(lCIF$COOR$VAL[1, "label"])
    ## print(lCIF$COOR$VAL[2, "label"])
  }

  ## The selection looks like this:
  ##
  ##         id index
  ## [1,] 12649    32
  ## [2,] 12649   131
  ## [3,] 12649   234
  ##
  ## > selection[1,][1]
  ## id
  ## 12649

  ## Set package global variables
  pkg$selection <- selection

  ## Restore the device ID.
  rgl::set3d(cur.dev, silent = TRUE)

  return(ret)
}
