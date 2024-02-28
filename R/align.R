## This file is part of the rgl.cry package
##
## Functions

#' Align crystal and diffraction pattern
#'
#' Align crystal and diffraction pattern and displayed.
#'
#' There is no z-axis alignment support because the visualization was created
#' with the analogy of selected area electron diffraction (SAED) on transmission
#' electron microscope (TEM) which typically have up to two axes.  However you
#' can rotate around the z-axis by the drag originates near the window edge.
#'
#' @param ax An axis to align
#' @param dev RGL device to apply.  Defaults to current device.
#' @param verbose logical: Should the report be suppressed?
#'
#' @return No return value, called for side effects.
#'
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'  align("a")
#'  align("rb")
#'  align("1 1 0")
#'  align("60 -30")
#'  align(dev = 123, "a")
#' }
#' }
align <- function(ax, dev = NULL, verbose = TRUE) {

  list(dev = dev, ax = ax, verbose = verbose)

  ## Select device
  if (missing(dev)) {
    tgt.dev <- rgl::cur3d()
  } else {
    tgt.dev <- dev
  }

  if (missing(ax)) {
    ax <- "empty"
  }


  ## Get the current pair of demo and check if it exist.
  inst <- pkg$inst
  idx <- which(inst$dp.dev == tgt.dev | inst$cry.dev == tgt.dev)
  ## Since integer(0) is returned when the element is empty, it will be judged
  ## by length().

  if (length(idx) == 0) {
    stop("The device was lost.\n")
  }


  ## Get the parameters.
  uc <- inst[[idx, "uc"]]
  ruc <- inst[[idx, "ruc"]]
  cry.dev <- inst[[idx, "cry.dev"]]
  cry.widget.id <- inst[[idx, "cry.widget.id"]]
  cry.root.id <- inst[[idx, "cry.root.id"]]
  cry.panel.id <- inst[[idx, "cry.panel.id"]]
  dp.dev <- inst[[idx, "dp.dev"]]
  dp.widget.id <- inst[[idx, "dp.widget.id"]]
  dp.root.id <- inst[[idx, "dp.root.id"]]
  dp.panel.id <- inst[[idx, "dp.panel.id"]]
  drawDp <- inst[[idx, "drawDp"]]
  ## drawCry <- ...

  ## Get the userMatrix of the beginning.
  cur.dev <- rgl::cur3d() # Save the current device.
  if (!is.na(cry.dev)) {
    rgl::set3d(cry.dev, silent = TRUE)
    rgl::useSubscene3d(cry.panel.id)
  } else if (!is.na(dp.dev)) {
    rgl::set3d(dp.dev, silent = TRUE)
    rgl::useSubscene3d(dp.panel.id)
  }
  umatPre <- rgl::par3d("userMatrix")
  rgl::set3d(cur.dev, silent = TRUE) # Restore the current device.


  ## Identity matrix.
  umatX <- rgl::identityMatrix()
  umatY <- rgl::identityMatrix()
  umatZ <- rgl::identityMatrix()


  ## ------------------------------------------------------------

  fromPrevisou <- FALSE

  if (length(grep("^[rabc]{1,2}$", ax)) == 1) { # a, b, c, ra, rb, rc

    ## Align the specified lattice vector to be perpedicular the screen (z).

    xyzf <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE) #

    if (length(grep("^r", ax)) != 1) { #
      xyz <- cry::frac_to_orth(
        xyzf, uc$a, uc$b, uc$c,
        uc$alpha, uc$beta, uc$gamma, 2
      )
      e <- matrix(unlist(xyz), ncol = 3, byrow = FALSE)
    } else {
      xyzr <- cry::frac_to_orth(
        xyzf, ruc$ar, ruc$br, ruc$cr,
        ruc$alpha, ruc$beta, ruc$gamma, 2
      )
      e <- matrix(unlist(xyzr), ncol = 3, byrow = FALSE)
    }

    ax <- gsub("r", "", ax)

    ## Rotate the target axis to face the line of sight.

    if (ax == "a") {
      t1 <- 1 # target axis
      t2 <- 2 # second target rotate around the target to align x.
    } else if (ax == "b") {
      t1 <- 2
      t2 <- 3
    } else if (ax == "c") {
      t1 <- 3
      t2 <- 1
    }

    tmp <- e[t1, ] * c(0, 1, 1) # Project the axis to the yz-plane.
    rotX <- acos(tmp[3] / norm(tmp, "2")) # Angle between projected tgt. and z.
    umatX <- rgl::rotationMatrix(rotX, 1, 0, 0) # rot. mat. around the x axis.
    e <- e %*% umatX[1:3, 1:3] # Update the vector's coordinates.

    tmp <- e[t1, ] * c(1, 0, 1) # Project the axis to the xz-plane
    rotY <- acos(tmp[3] / norm(tmp, "2"))
    umatY <- rgl::rotationMatrix(rotY, 0, 1, 0) # rotation around the y axis
    e <- e %*% umatY[1:3, 1:3] # apply rotate

    tmp <- e[t2, ] * c(1, 1, 0) # Project the axis to the xy-plane
    rotZ <- acos(tmp[1] / norm(tmp, "2"))
    rotZ <- ifelse(tmp[2] > 0, rotZ, -rotZ) #
    umatZ <- rgl::rotationMatrix(rotZ, 0, 0, 1)
  } else if (length(grep("^[-]?[0-9]+ [-]?[0-9]+ [-]?[0-9]+$", ax)) == 1) {
    ## Align the specified lattice plane (hkl) parallel to the screen.

    xyzf <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE) #
    xyzf <- do.call(as.numeric, strsplit(ax, " "))
    xyzr <- cry::frac_to_orth(
      xyzf, ruc$ar, ruc$br, ruc$cr,
      ruc$alpha, ruc$beta, ruc$gamma, 2
    )
    e <- matrix(unlist(xyzr), ncol = 3, byrow = FALSE)

    tmp <- e * c(0, 1, 1) # Project the axis to the yz-plane.
    rotX <- acos(tmp[3] / norm(tmp, "2")) # Angle between projected tgt. and z.
    umatX <- rgl::rotationMatrix(rotX, 1, 0, 0) # rot. mat. around the x axis.
    e <- e %*% umatX[1:3, 1:3] # Update the vector's coordinates.

    tmp <- e * c(1, 0, 1) # Project the axis to the xz-plane
    rotY <- acos(tmp[3] / norm(tmp, "2"))
    umatY <- rgl::rotationMatrix(rotY, 0, 1, 0) # rotation around the y axis
    e <- e %*% umatY[1:3, 1:3] # apply rotate
  } else if (length(grep("^[-]?[0-9.]+ [-]?[0-9.]+$", ax)) == 1) { # rotX, rotY

    ## Perform the specified X and Y rotation.
    ## - (Rotation from the initial rotation matrix state.)
    ## + (Rotation from the current rotation matrix state.)

    ax <- gsub("^[ ]+", "", ax)
    ax <- gsub("[ ]+", " ", ax)

    rot <- do.call(as.numeric, strsplit(ax, " "))
    rotX <- rot[1] * pi / 180
    rotY <- rot[2] * pi / 180

    umatX <- rgl::rotationMatrix(rotX, 1, 0, 0) # rot. mat. around the x axis.
    umatY <- rgl::rotationMatrix(rotY, 0, 1, 0) # rot. around the y axis
    umatZ <- rgl::rotationMatrix(0, 0, 0, 1) # rot. around the z axis

    fromPrevisou <- TRUE
  } else {
    stop("Syntax error.\n")
  }


  ## The rotation matrix depends on the state to which it is applied.
  if (fromPrevisou == TRUE) {
    umat <- solve(umatZ) %*% solve(umatY) %*% solve(umatX) %*% umatPre
  } else {
    umat <- solve(umatZ) %*% solve(umatY) %*% solve(umatX)
  }



  ## ------------------------------------------------------------
  ## Apply userMatrix to the scenes.

  ## Save the current device.
  cur.dev <- rgl::cur3d()

  if (!is.na(dp.dev)) {
    rgl::set3d(dp.dev, silent = TRUE)
    scenes <- list(dp.panel.id, dp.root.id, dp.widget.id)

    lapply(scenes, function(v) {
      rgl::par3d(subscene = v, userMatrix = umat)
    })

    rgl::useSubscene3d(dp.panel.id)
    drawDp()
  }

  if (!is.na(cry.dev)) {
    rgl::set3d(cry.dev, silent = TRUE)
    scenes <- list(cry.panel.id, cry.root.id, cry.widget.id)

    lapply(scenes, function(v) {
      rgl::par3d(subscene = v, userMatrix = umat)
    })

    ## rgl::useSubscene3d(dp.panel.id)
    ## drawDp()
  }

  ## Restore the current device.
  rgl::set3d(cur.dev, silent = TRUE)



  ## ------------------------------------------------------------

  ## Print out the diff to previous.
  ##
  ## Memo
  ##
  ##   Rx <- matrix(c( 1,  0,  0,  0, Cx, -Sx,   0, Sx, Cx), ncol=3, byrow=T)
  ##   Ry <- matrix(c(Cy,  0, Sy,  0,  1,   0, -Sy,  0, Cy), ncol=3, byrow=T)
  ##   Rz <- matrix(c(Cz, -Sz, 0, Sz, Cz,   0,   0,  0,  1), ncol=3, byrow=T)
  ##
  ##   R <- Rz Ry Rx
  ##
  ##        CyCz, CzSx-CxSz, CxCzSy+SySz,
  ##        CySz, SxSz+CxCz, CxSySz-CzSx,
  ##         -Sy,      CySx,        CxCy
  ##
  ##   θy = asin(R[3,1])         from asin(-Sy)
  ##   θx = atan(R[3,2]/R[3,3])  from atan(CySx/CxCy)  from atan(Sx/Cx)
  ##   θz = atan(R[2,1]/R[1,1])  from atan(CySz/CyCz)  from atan(Sz/Cz)
  ##

  if (verbose == TRUE) {

    ## diff <- solve(umatPre) %*% umat
    diff <- umat %*% solve(umatPre)
    rotX <- -atan(diff[3, 2] / diff[3, 3]) * 180 / pi
    rotY <- -asin(-diff[3, 1]) * 180 / pi
    rotZ <- -atan(diff[2, 1] / diff[1, 1]) * 180 / pi

    ## To remove sign such as -0.  This is consistent with the %.0f later.
    rotX <- ifelse(abs(rotX) < 1, 0, rotX)
    rotY <- ifelse(abs(rotY) < 1, 0, rotY)
    rotZ <- ifelse(abs(rotZ) < 1, 0, rotZ)

    ## This is the result of the following rotations from previous.
    message(sprintf(
      "Current state from previous: x y z (deg): %.0f %.0f %.0f (deg)",
      rotX, rotY, rotZ
    ))
  }

}
