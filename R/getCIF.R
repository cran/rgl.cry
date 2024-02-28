## This file is part of the rgl.cry package
##
## Functions

#' Examples of using the cry and rgl packages together.
#'
#' Output a lCIF.
#' 
#' @param dev The device that is used to extract the lCIF.  The default is
#' current device.
#'
#' @return A named list, the same as that of 'cry::readCIF.'
#'
#' @export
#' @examples
#' \donttest{
#' getCIF()
#' }
getCIF <- function(dev = NULL) {
  list(dev = dev)

  if (is.null(dev)) {
    cur.dev <- rgl::cur3d()
  } else {
    cur.dev <- dev
  }

  inst <- pkg$inst # Get the current list of instance.
  idx <- which(inst$dp.dev == cur.dev | inst$cry.dev == cur.dev)
  ## Since integer(0) is returned when the element is empty, it will be judged
  ## by length().

  if (length(idx) == 0) {
    ret <- NULL
  } else {
    ret <- inst[[idx, "lCIF"]]
  }

  return(ret)
}
