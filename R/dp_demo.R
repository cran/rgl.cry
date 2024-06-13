## This file is part of the rgl.cry package
##
## Functions and variables:

#' Examples of using the cry and rgl packages together.
#'
#' Read a file in CIF formats, set and the parameters, calculates them, draws
#' the reciprocal lattice map with a cell widget.
#'
#' If no file argument is provided, and `cry_demo()` has been opened without
#' paired `dp_demo()`, the CIF parameters of already opened `cry_demo()` will be
#' used.
#'
#' Interactive rotation, zooming, and panning of structures are possible using
#' the 3D graphics library \code{rgl}.  When the drag originates near the window
#' edge (within 5%), perform a Z-axis rotation.
#'
#' This function also performs powder diffraction simulation and saves the
#' results to a file in the working directory.  Currently, it doesn't account
#' for atomic ionization and uses standard atomic scattering factors.
#'
#' @param file Optional file in CIF formats.
#' @param reso A real number. The highest data resolution, in angstroms.  If the
#' default value takes a long time to process displaying due to the large number
#' of lattice points, you can expect to improve performance by increasing the
#' value.
#' @param ews.r Ewald sphere radius in angstrom^-1.
#' @param zoom A positive value indicating the current scene magnification.
#' @param xrd A logical value indicating whether to create an X-ray diffraction
#' pattern simulation result file.
#'
#' @return An integer the device number of the currently window.
#'
#' @export
#'
#' @examples
#' dp_demo()
#' dp_demo(system.file("orthorhombic_p.cif", package = "rgl.cry"))
#' dp_demo(system.file("orthorhombic_p.cif", package = "rgl.cry"), res = 2.0)
#'
#' \donttest{
#' if (interactive()) {
#'  dp_demo(file, zoom = 0.5)
#'  dp_demo("https://www.crystallography.net/cod/foo.cif")
#' }
#' }
dp_demo <- function(file = NULL, reso = 1.2, ews.r = 40, zoom = 0.5, xrd = FALSE) {

  list(file = file, reso = reso, ews.r = ews.r, zoom = zoom, xrd = xrd)

  ## File or lCIF object to use.
  if (!is.null(file)) {
    lCIF <- cry::readCIF(file)
  } else {
    lCIF <- getCIF()
    if (is.null(lCIF)) {
      file <- system.file("monoclinic_a.cif", package = "rgl.cry")
      lCIF <- cry::readCIF(file)
    }
  }


  ## Extract the crystal parameters from a lCIF.
  a <- lCIF$HEADER$CELL$A$VAL
  b <- lCIF$HEADER$CELL$B$VAL
  c <- lCIF$HEADER$CELL$C$VAL
  aa <- lCIF$HEADER$CELL$ALPHA$VAL # degree
  bb <- lCIF$HEADER$CELL$BETA$VAL
  cc <- lCIF$HEADER$CELL$GAMMA$VAL
  SG <- cry::findHM(lCIF$HEADER$HM)

  uc <- cry::create_unit_cell(a, b, c, aa, bb, cc)
  ruc <- cry::create_rec_unit_cell(uc)

  ## Create a data of Miller indices.
  ## 
  ## The systematic absences have been taken into account in
  ## cry::generate_miller() by calling cry::deplete_systematic_absences().
  hkl <- cry::generate_miller(uc, SG, reso) # list the h+k+l <= reso.

  xyzf <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE) #
  xyz <- cry::frac_to_orth(
    xyzf, a, b, c,
    aa, bb, cc, 2
  )
  ea1 <- as.numeric(xyz[1, ])
  ea2 <- as.numeric(xyz[2, ])
  ea3 <- as.numeric(xyz[3, ])

  V <- as.numeric(ea1 %*% pracma::cross(ea2, ea3))
  eb1 <- pracma::cross(ea2, ea3) / V
  eb2 <- pracma::cross(ea3, ea1) / V
  eb3 <- pracma::cross(ea1, ea2) / V

  pos <- t(apply(hkl, 1, function(v) {
    v["H"] * eb1 + v["K"] * eb2 + v["L"] * eb3
  }))

  ## Ewald sphere and text label settings.
  ##ews.r <- 40 # relativistic wave number of electron beam at 200 kV
  ews.pos <- c(0, 0, ews.r) # center of Ewald sphere.
  text.offset <- c(0.03, 0.03, 0) # offset for text.



  ## ------------------------------------------------------------
  ## Prepare the coordinates of atoms in the unit cell.

  num <- dim(lCIF$COOR$VAL)[1] # Number of non-equivalent atoms.
  sites <- NULL # The label, coordinate and name of that atom.
  for (i in seq_len(num)) {
    sites <- rbind(
      sites,
      c(
        lCIF$COOR$VAL[i, "label"],
        lCIF$COOR$VAL[i, "fract_x"],
        lCIF$COOR$VAL[i, "fract_y"],
        lCIF$COOR$VAL[i, "fract_z"],
        gsub("[[:digit:][:punct:]]+", "", lCIF$COOR$VAL[i, "label"]) # name
      )
    )
  }

  ## Coordinates of equivalent points for this space group.
  ##
  ##   > pos.xyz
  ##        [,1] [,2] [,3]
  ##   [1,] "x"  "y"  "z"
  ##   [2,] "-x" "-y" "z"
  ##   [3,] "-x" "-y" "-z"
  ##   [4,] "x"  "y"  "-z"
  ##
  op_xyz_list <- cry::syminfo_to_op_xyz_list(cry::findHM(lCIF$HEADER$HM))
  xyz <- op_xyz_list$SYMOP # get _symmetry_equiv_pos_as_xyz
  pos.xyz <- NULL
  for (i in seq_len(length(xyz))) {
    ## Split the string by commas to obtain individual parts then unlist.
    pos.xyz <- rbind(pos.xyz, unlist(strsplit(xyz[i], ","))) #
  }

  ## Assign the atomic coordinates to the equivalent points for each atom.
  ##
  ##   > pos.cry
  ##   [[1]]
  ##        [,1] [,2] [,3]
  ##   [1,]  0.5    0    0
  ##   [2,] -0.5    0    0
  ##   [3,] -0.5    0    0
  ##   [4,]  0.5    0    0
  ##
  ##   [[2]]
  ##           [,1]    [,2] [,3]
  ##   [1,]  0.1645  0.3363    0
  ##   [2,] -0.1645 -0.3363    0
  ##   [3,] -0.1645 -0.3363    0
  ##   [4,]  0.1645  0.3363    0
  ##
  pos.cry <- NULL
  for (i in seq_len(num)) {
    var <- list(
      x = as.numeric(sites[i, 2]),
      y = as.numeric(sites[i, 3]),
      z = as.numeric(sites[i, 4])
    )
    pos.cry <- c(pos.cry, list(apply(
      pos.xyz, c(1, 2),
      function(v) eval(parse(text = v), var)
    )))
  }

  ## Add the offset to each lattice point.
  ##
  ##   > xyz2
  ##        [,1] [,2]    [,3]
  ##   [1,] "x"  "y"     "z"
  ##   [2,] "x"  "y+1/2" "z+1/2"
  ##   > sft
  ##        [,1] [,2] [,3]
  ##   [1,]    0  0.5  0.5
  ##
  cenop <- op_xyz_list$CENOP

  xyz2 <- NULL
  sft <- NULL
  for (i in seq_len(length(cenop))) {
    xyz2 <- rbind(xyz2, unlist(strsplit(cenop[i], ",")))
  }
  sft <- apply(
    xyz2, c(1, 2),
    function(v) eval(parse(text = v), list(x = 0, y = 0, z = 0))
  )
  sft <- sft[apply(sft, 1, function(x) any(x != 0)), ] # remove 0, 0, 0
  ## If sft has only one row, the output is a vector, so a workaround is
  ## needed.
  if (length(sft) != 0) {
    sft <- t(as.matrix(sft))
    for (i in seq_len(num)) {
      tmp <- matrix(apply(pos.cry[[i]], 1, function(v) {
        v + sft
      }), ncol = 3, byrow = TRUE)
      pos.cry[[i]] <- rbind(pos.cry[[i]], tmp)
    }
  }

  ## Add 1 to the atomic coordinate if it is negative to make it within the
  ## unit cell.
  for (i in seq_len(num)) {
    pos.cry[[i]] <- t(apply(pos.cry[[i]], 1, function(v) {
      unlist(lapply(v, function(w) ifelse(w < 0, w + 1, w)))
    }))
  }

  ## Similarly, if the value is greater than 1, subtract 1.
  for (i in seq_len(num)) {
    pos.cry[[i]] <- t(apply(pos.cry[[i]], 1, function(v) {
      unlist(lapply(v, function(w) ifelse(w > 1, w - 1, w)))
    }))
  }

  ## Remove the same coordinates.
  for (i in seq_len(num)) {
    pos.cry[[i]] <- unique(pos.cry[[i]])
  }



  ## ------------------------------------------------------------
  ## Calculate the structure factor.
  ##
  ## This is also used for the powder X-ray diffraction pattern simulation.

  ## The intensity of powder X-ray diffraction is dominated by the structure
  ## factor |F|^2 as the following equation [1, Chap. 16]:
  ##
  ##   F(hkl) = Σj f_j exp(2πi(h x_j + k y_j + l z_j)),
  ##
  ## where f_j is the atomic scattering factor [1, Chap. 3.7] and I'm refering
  ## [2] for the caliculation but heavily modified when using.  In perticular
  ## the value is fitted using smooth.spline.  Also we need modification when it
  ## is applyed to the electron diffraction [1, eq. 3.9].
  ##
  ## Here, the scattering factor is summed over all atoms in the Bravais
  ## lattice, so the multiplicity (often seen in explanations) is not used in
  ## this R script.
  ##
  ## The Lorentz–polarization correction factor Lp is multiplied as a correction
  ## term for the intensity.  The Lorentz-polarization correction factor is as
  ## follows [4]:
  ##
  ##   (1 + A (cos 2θ)^2) / ((1 + A)(sin 2θ)),
  ##
  ## where A = (cos 2θm)^2 and θm is the Bragg angle of the monochromator
  ## crystal.  However, since A is unknown at this point, we set it to 1.  In
  ## terms of Lorentz correction, we use the conventional form of
  ## 1/((sin θ)^2 (cos θ)) instead of 1/sin 2θ, we now obtain the following
  ## equation:
  ##
  ##   (1 + (cos 2θ)^2) / ((sin θ)^2 (cos θ))
  ##

  ## Atomic scattering factor.

  ScatFac <- NULL
  load(system.file("extdata/sc.Rda", package = "rgl.cry"))

  ## Lorentz–polarization correction factor.
  lp <- function(th) {
    (1 + cos(2 * th)^2) / (sin(th)^2 * cos(th))
  }

  StrFac <- NULL # structure factor

  for (j in seq_len(dim(hkl)[1])) {
    h <- hkl[j, ]$H
    k <- hkl[j, ]$K
    l <- hkl[j, ]$L

    f <- 0
    for (i in seq_len(num)) {
      ## The label is used to select the atom.  We don't check the ionization of
      ## atom and ignore now.
      element <- gsub("[[:digit:][:punct:]]+", "", lCIF$COOR$VAL[i, "label"])

      fx <- stats::predict(ScatFac[[element]], x = hkl[j, ]$S / 2)$y

      f <- f + fx * sum(apply(pos.cry[[i]], 1, function(v) {
        x <- v[1]
        y <- v[2]
        z <- v[3]
        complex(0, cos(2 * pi * (h * x + k * y + l * z)),
                sin(2 * pi * (h * x + k * y + l * z)))
      }))
    }

    f <- abs(f)
    StrFac[j] <- ifelse(f < 1e-6, 0, f)
  }

  pdp <- data.frame(
    h = hkl$H,
    k = hkl$K,
    l = hkl$L,
    d = 1 / hkl$S,
    absF = abs(StrFac),
    lp = lp(asin(0.7709167 * hkl$S)), # Assuming the use of Cu Kα 1.541833
    twotheta = asin(0.7709167 * hkl$S) * 360 / pi
  )

  ## write the data.
  ## with this data, we can plot the calculated powder X-ray diffraction
  ## pattern.

  if (xrd) {
    xrdfile <- paste0(
      "rgl.cry.dp.demo.",
      format(Sys.time(), "%Y-%m-%d_%H%M%S.dat")
    )
    sink(xrdfile, append = FALSE)
    print(pdp)
    sink()
  }


  ## apply the absF for the reciprocal lattice point radius.




  ## ------------------------------------------------------------
  ## Define a drawing function.

  drawDp <- function() {
    umat <- rgl::par3d("userMatrix")

    ## Remove previous plots (convenience objects, reciprocal lattice, label)
    rgl::pop3d(tag = c("rlpoints0", "rlpoints1", "rlpoints2"))

    ## Caliculate the lattice points that are within a certain distance of the
    ## Ewald sphere surface.

    ews.pos.new <- ews.pos %*% umat[1:3, 1:3]

    dist <- apply(pos, 1, function(v) {
      round(dist(rbind(ews.pos.new, v)), digits = 6)
    })
    vis <- lapply(dist, function(w) {
      ## Scale the range from 0-0.0333 to 0-1 then subtract from 1.
      v <- 1 - 30 * abs(w - ews.r)
      v <- ifelse(v < 0, 0, v) # Flooring the result to 0.
    })

    ## Place the sphere to prevent the draw area from being clipped.
    rgl::spheres3d(frame, r = 0, color = "green", alpha = 0, tag = "rlpoints0")

    ## Draw the reciprocal lattice points that satisfy the criteria, and set the
    ## alpha value of the lattice points that do not satisfy the criteria to 0.
    rgl::spheres3d(pos,
      r = 0.01, color = "black",
      alpha = vis, tag = "rlpoints1"
    )
    str <- paste(hkl$H, hkl$K, hkl$L)
    rgl::text3d(sweep(pos, 2, (text.offset %*% umat[1:3, 1:3]), "+"),
      texts = str, cex = 0.8, col = "blue", alpha = vis, tag = "rlpoints2"
    )

    ##
    ## The text label for hkl index placed each reciprocal lattice points with
    ## offset. To ensure that the text labels are always facing the viewer, the
    ## coordinate of label is counter-rotated around the bbox center, and then
    ## translated to each reciprocal points.
    ##
    ## Note:  At first, I made a mistake by doing the following, so the label
    ##        rotated slightly with rotation. The following does not add a
    ##        vector to each rows.
    ##
    ##        as.numeric(text.offset %*% umat[1:3,1:3]) + pos
    ##
  }



  ## ------------------------------------------------------------

  ## Open device.
  dp.dev <- rgl::open3d()

  inst <- pkg$inst
  if (is.na(inst[nrow(inst), "dp.dev"])) {
    inst[nrow(inst), "dp.dev"] <- dp.dev
  } else {
    inst[nrow(inst) + 1, "dp.dev"] <- dp.dev
  }

  ## I intended to share the behavior of the mouse drag operation with the
  ## cry_demo and dp_demo pairs that handle the same lCIF.  It's good that a
  ## pair will be made upon function invocation, but I haven't yet finalize the
  ## implementation details for that behavior.


  ## Delete child scence and reset the callback to the default.
  try(
    {
      while (length(rgl::subsceneInfo()$parent) != 0) {
        rgl::useSubscene3d(rgl::subsceneInfo()$parent)
      }
      if (length(rgl::subsceneInfo()$children) != 0) {
        rgl::delFromSubscene3d(ids = rgl::subsceneInfo()$children)
      }
    },
    silent = TRUE
  )


  ## initial view settings
  rgl::par3d(mouseMode = rgl::r3dDefaults$mouseMode)
  rgl::clear3d()
  rgl::par3d(windowRect = c(0, 0, 500, 500))
  rgl::view3d(theta = 0, phi = 0, zoom = 1, fov = 0)

  frame <- rbind(
    c(-1, -1, -1), c(-1, -1, 1), c(1, -1, -1), c(-1, 1, -1),
    c(-1, 1, 1), c(1, -1, 1), c(1, 1, -1), c(1, 1, 1)
  )


  ## ------------------------------------------------------------
  ## Subscene for main drawing.

  ## Save the scene id and set the mouseMode for this scene.
  dp.root.id <- rgl::currentSubscene3d() # subsceneInfo()$id
  rgl::par3d(mouseMode = c("none", "none", "none", "none", "none"))


  ## ------------------------------------------------------------
  ## Subscene for widget

  rgl::newSubscene3d(newviewport = c(0, 0, 150, 150), mouseMode = "replace")

  ## Save the scene id and set the mouseMode for this scene.
  dp.widget.id <- rgl::currentSubscene3d() # subsceneInfo()$id
  rgl::par3d(mouseMode = c("none", "none", "none", "none", "none"))

  ## Place the dummy sphere to prevent the draw area from modification.
  ## r=0 is prevention of protrusion.
  rgl::spheres3d(frame, r = 0, color = "green", alpha = 0) #

  ## center
  oo <- c(0, 0, 0)
  rgl::spheres3d(oo, r = 0.05, color = "black") #


  ## lattice basis vectors
  xyzf <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE) #


  ## lattice
  ## 
  ## When ochoice = 2, as the help document of frac_to_orth() states, "The X axis
  ## is along a*; the Y axis lies in the (a*,b*) plane; the Z axis is,
  ## consequently, along c."
  ## 
  xyz <- cry::frac_to_orth(
    xyzf, uc$a, uc$b, uc$c,
    uc$alpha, uc$beta, uc$gamma, 2
  )
  xyz_max <- max(xyz[1, ], xyz[2, ], xyz[3, ])
  ea1 <- as.numeric(xyz[1, ]) / xyz_max
  ea2 <- as.numeric(xyz[2, ]) / xyz_max
  ea3 <- as.numeric(xyz[3, ]) / xyz_max

  lines <- rbind(
    oo, ea1, oo, ea2, ea1 + ea2, ea1, ea1 + ea2, ea2,
    ea3, ea1 + ea3, ea3, ea2 + ea3, ea1 + ea2 + ea3, ea1 + ea3, ea1 + ea2 + ea3, ea2 + ea3,
    oo, ea3, ea1, ea1 + ea3, ea2, ea2 + ea3, ea1 + ea2, ea1 + ea2 + ea3
  )
  rgl::segments3d(lines, col = "gray", lwd = 1.5)
  rgl::text3d(ea1 * 1.2, texts = "a", cex = 1, col = "gray")
  rgl::text3d(ea2 * 1.2, texts = "b", cex = 1, col = "gray")
  rgl::text3d(ea3 * 1.2, texts = "c", cex = 1, col = "gray")


  ## reciprocal lattice
  V <- as.numeric(ea1 %*% pracma::cross(ea2, ea3))
  eb1 <- pracma::cross(ea2, ea3) / V
  eb2 <- pracma::cross(ea3, ea1) / V
  eb3 <- pracma::cross(ea1, ea2) / V
  eb_max <- max(eb1, eb2, eb3)
  eb1 <- eb1 / eb_max
  eb2 <- eb2 / eb_max
  eb3 <- eb3 / eb_max

  lines <- rbind(
    oo, eb1, oo, eb2, eb1 + eb2, eb1, eb1 + eb2, eb2,
    eb3, eb1 + eb3, eb3, eb2 + eb3, eb1 + eb2 + eb3, eb1 + eb3, eb1 + eb2 + eb3, eb2 + eb3,
    oo, eb3, eb1, eb1 + eb3, eb2, eb2 + eb3, eb1 + eb2, eb1 + eb2 + eb3
  )
  rgl::segments3d(lines, col = "lightpink", lwd = 1.5)
  rgl::text3d(eb1 * 1.1, texts = "a*", cex = 1, col = "lightpink")
  rgl::text3d(eb2 * 1.1, texts = "b*", cex = 1, col = "lightpink")
  rgl::text3d(eb3 * 1.1, texts = "c*", cex = 1, col = "lightpink")


  ## ------------------------------------------------------------
  ## Subscene on the top of scenes for event handling.

  rgl::newSubscene3d(newviewport = c(0, 0, 500, 500), mouseMode = "replace")

  ## Save the scene id and set the mouseMode for this scene.
  dp.panel.id <- rgl::currentSubscene3d()
  rgl::par3d(mouseMode = c("none", "user", "none", "none", "pull"))



  ## ------------------------------------------------------------
  ## Mouse event handling.

  start <- list()
  start$time <- 0

  begin <- function(x, y) {

    ## Save parameters.
    time.current <- as.numeric(Sys.time()) * 1000 # micro to milli
    time.difference <- time.current - start$time
    start$time <<- time.current

    ## These parameters were saved when this function was set as a callback and
    ## are not chaged when another function call performed:
    ##   dp.dev, dp.widget.id, dp.root.id and dp.panel.id

    rgl::set3d(dp.dev, silent = TRUE)
    rgl::useSubscene3d(dp.panel.id)
    umat <- rgl::par3d("userMatrix")


    ## Retrieve the cry and dp pair of this instance.
    inst <- pkg$inst # Get the current list of instance.
    idx <- which(inst$dp.dev == dp.dev)
    start$cry.dev <<- inst[inst$dp.dev == dp.dev, "cry.dev"]
    start$cry.widget.id <<- inst[idx, "cry.widget.id"]
    start$cry.root.id <<- inst[idx, "cry.root.id"]
    start$cry.panel.id <<- inst[idx, "cry.panel.id"]


    ## The rotation is reset to its original value when the mouse is
    ## double-clicked.
    if (time.difference > 100 && time.difference < 300) {
      rgl::pop3d(tag = c("rlpoints0", "rlpoints1", "rlpoints2"))
      umat <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
      rgl::par3d(subscene = dp.widget.id, userMatrix = umat)
      rgl::par3d(subscene = dp.root.id, userMatrix = umat)
      rgl::par3d(subscene = dp.panel.id, userMatrix = umat)
      drawDp()

      ## Update cry_demo drawings if exists.
      if (!is.na(start$cry.dev)) {
        rgl::set3d(start$cry.dev, silent = TRUE)
        rgl::par3d(subscene = start$cry.widget.id, userMatrix = umat)
        rgl::par3d(subscene = start$cry.root.id, userMatrix = umat)
        rgl::par3d(subscene = start$cry.panel.id, userMatrix = umat)
        rgl::useSubscene3d(start$cry.panel.id)
        ## drawCry() # not necessary.
        rgl::set3d(dp.dev, silent = TRUE)
      }
    }

    ## Save the current mouse cursor position.
    start$x <<- x
    start$y <<- y
    start$umat <<- umat
  }

  ## Rotate and redraw by draging the mouse.
  update <- function(x, y) {

    ## Get the umat at the begining.
    umat <- start$umat # call begin then call update without sequencially
    viewport <- rgl::par3d("viewport", dev = dp.dev)

    w <- viewport[["width"]]
    h <- viewport[["height"]]
    x <- (x - start$x) / w / 1 # 1 is empirically derived.
    y <- (y - start$y) / h / 1

    rot <- 0
    if (start$x > 0.95 * w) {
      rot <- -y
    } else if (start$x < 0.05 * w) {
      rot <- y
    } else if (start$y > 0.95 * h) {
      rot <- x
    } else if (start$y < 0.05 * h) {
      rot <- -x
    }

    if (rot != 0) {
      ## When the drag originates near the window edge (within 5%),
      ## perform a Z-axis rotation.
      ez <- (c(0, 0, 1, 1) %*% umat)[1:3]
      umat <- umat %*% rgl::rotationMatrix(rot, ez[1], ez[2], ez[3])
    } else {
      ## When the drag on the window, rotate the axis that is perpendicular to
      ## the drag direction.
      ex <- (c(1, 0, 0, 1) %*% umat)[1:3]
      ey <- (c(0, 1, 0, 1) %*% umat)[1:3]
      umat <- umat %*% rgl::rotationMatrix(x, ey[1], ey[2], ey[3])
      umat <- umat %*% rgl::rotationMatrix(y, ex[1], ex[2], ex[3])
    }

    ## Control the scene using this instance's device ID.
    rgl::set3d(dp.dev, silent = TRUE)
    rgl::par3d(subscene = dp.widget.id, userMatrix = umat)
    rgl::par3d(subscene = dp.root.id, userMatrix = umat)
    rgl::par3d(subscene = dp.panel.id, userMatrix = umat)
    rgl::useSubscene3d(dp.panel.id)
    drawDp()

    ## If a pair of cry_demo exists, perform the same operation.
    if (!is.na(start$cry.dev)) {
      rgl::set3d(start$cry.dev, silent = TRUE)
      rgl::par3d(subscene = start$cry.widget.id, userMatrix = umat)
      rgl::par3d(subscene = start$cry.root.id, userMatrix = umat)
      rgl::par3d(subscene = start$cry.panel.id, userMatrix = umat)
      rgl::useSubscene3d(start$cry.panel.id)
      ## drawCry() # not necessary.
      rgl::set3d(dp.dev, silent = TRUE)
    }

  }

  rotate <- function(r) {
    ## Currenly not used.
    ## zoom <- rgl::par3d()$zoom + ifelse(r == 1, 0.015, -0.015)
    ## rgl::par3d(subscene = dp.widget.id, zoom = zoom)
    ## rgl::par3d(subscene = dp.root.id, zoom = zoom)
    ## rgl::par3d(subscene = dp.panel.id, zoom = zoom)
    ## print(zoom)
  }



  ## ------------------------------------------------------------
  ## Start.

  ## Set the callback for the scene dp.panel.id.
  rgl::useSubscene3d(dp.panel.id)
  rgl::rgl.setMouseCallbacks(1, begin, update)
  ## rgl.setWheelCallback(rotate)

  ## Initial view setting.
  rgl::par3d(zoom = zoom)

  ##
  ## cat("Use this panel ID to change the view settings: ", dp.panel.id)

  ##
  drawDp()

  ## Set package global variables
  inst[[nrow(inst), "uc"]] <- uc
  inst[[nrow(inst), "ruc"]] <- ruc
  inst[[nrow(inst), "hkl"]] <- hkl
  inst[[nrow(inst), "lCIF"]] <- lCIF
  inst[[nrow(inst), "drawDp"]] <- drawDp
  inst[nrow(inst), "dp.root.id"] <- dp.root.id
  inst[nrow(inst), "dp.widget.id"] <- dp.widget.id
  inst[nrow(inst), "dp.panel.id"] <- dp.panel.id
  assign("inst", inst, pkg)
  assign("drawDp", drawDp, pkg)

  ##
  return(as.numeric(dp.dev))
}



## ------------------------------------------------------------
## Reference
##
## [1] David B. Williams and C. Barry Carter, Transmission Electron
##     Microscopy, Springer US, 2009. doi: 10.1007/978-0-387-76501-3
##
## [2] T. Hanashima, "Contents of atomic scattering factors' table in S. Sasaki
##     (1987), KEK Report 87-3", constructed in web page in 2001.
##     https://www2.kek.jp/imss/pf/tools/sasaki/sinram/sinram.html (accessed
##     2024-02)
##
## [3] Relativistic Scattering Factors for Atoms and Atomic Valence Shells
##     http://harker.chem.buffalo.edu/group/ptable.html (accessed 2024-02)
##
## [4] International Union of Crystallography. "Lorentz–polarization
##     correction." In Online Dictionary of Crystallography. Accessed February
##     15, 2024.
##     https://dictionary.iucr.org/Lorentz%E2%80%93polarization_correction
##     (accessed 2024-02-14)
##


#' @noRd
.makeFittedAtomicScatteringFactor <- function(save = FALSE) {

  list(save = save)

  ScatFac <- list()
  files <- list("H~0.txt", "HE~0.txt", "LI~0.txt", "LI~p1.txt", "BE~0.txt", "BE~p2.txt", "B~0.txt", "B~p2.txt", "B~p3.txt", "C~0.txt", "N~0.txt", "O~0.txt", "O~m1.txt", "O~m2.txt", "O~m1.9.txt", "O~m1.4.txt", "F~0.txt", "F~m1.txt", "NE~0.txt", "NA~0.txt", "NA~p1.txt", "MG~0.txt", "MG~p1.txt", "MG~p2.txt", "MG~p1.9.txt", "AL~0.txt", "AL~p1.txt", "AL~p2.txt", "AL~p3.txt", "SI~0.txt", "SI~p1.txt", "SI~p2.txt", "SI~p3.txt", "SI~p4.txt", "P~0.txt", "S~0.txt", "S~m1.txt", "CL~0.txt", "CL~m1.txt", "AR~0.txt", "K~0.txt", "K~p1.txt", "CA~0.txt", "CA~p1.txt", "CA~p2.txt", "SC~0.txt", "SC~p1.txt", "SC~p3.txt", "TI~0.txt", "TI~p2.txt", "TI~p3.txt", "TI~p4.txt", "V~0.txt", "V~p2.txt", "V~p3.txt", "V~p5.txt", "CR~0.txt", "CR~p2.txt", "CR~p3.txt", "MN~0.txt", "MN~p1.txt", "MN~p2.txt", "MN~p3.txt", "MN~p4.txt", "FE~0.txt", "FE~p1.txt", "FE~p2.txt", "FE~p3.txt", "CO~0.txt", "CO~p1.txt", "CO~p2.txt", "CO~p3.txt", "CO~p1.2.txt", "NI~0.txt", "NI~p1.txt", "NI~p2.txt", "NI~p3.txt", "NI~p0.7.txt", "CU~0.txt", "CU~p1.txt", "CU~p2.txt", "ZN~0.txt", "ZN~p2.txt", "GA~0.txt", "GA~p3.txt", "GE~0.txt", "GE~p4.txt", "AS~0.txt", "SE~0.txt", "BR~0.txt", "BR~m1.txt", "KR~0.txt", "RB~0.txt", "RB~p1.txt", "SR~0.txt", "SR~p2.txt", "Y~0.txt", "Y~p3.txt", "ZR~0.txt", "ZR~p4.txt", "NB~0.txt", "NB~p3.txt", "NB~p5.txt", "MO~0.txt", "MO~p3.txt", "MO~p5.txt", "TC~0.txt", "RU~0.txt", "RH~0.txt", "PD~0.txt", "AG~0.txt", "AG~p1.txt", "CD~0.txt", "CD~p2.txt", "IN~0.txt", "IN~p3.txt", "SN~0.txt", "SN~p2.txt", "SB~0.txt", "SB~p3.txt", "TE~0.txt", "I~0.txt", "I~m1.txt", "XE~0.txt", "CS~0.txt", "CS~p1.txt", "BA~0.txt", "BA~p2.txt", "LA~0.txt", "CE~0.txt", "PR~0.txt", "ND~0.txt", "PM~0.txt", "SM~0.txt", "EU~0.txt", "GD~0.txt", "GD~p3.txt", "TB~0.txt", "DY~0.txt", "HO~0.txt", "ER~0.txt", "ER~p3.txt", "TM~0.txt", "YB~0.txt", "LU~0.txt", "HF~0.txt", "TA~0.txt", "W~0.txt", "W~p6.txt", "RE~0.txt", "OS~0.txt", "IR~0.txt", "PT~0.txt", "PT~p2.txt", "PT~p4.txt", "AU~0.txt", "AU~p1.txt", "AU~p3.txt", "HG~0.txt", "HG~p1.txt", "HG~p2.txt", "TL~0.txt", "PB~0.txt", "PB~p2.txt", "BI~0.txt", "BI~p3.txt", "U~0.txt")

  for (i in files) {
    file <- paste0("../scattering_factors/data/", i)

    df <- utils::read.table(
      file = file,
      header = FALSE,
      sep = "\t",
      fileEncoding = "UTF-16LE",
      comment.char = "#",
      blank.lines.skip = TRUE,
      skipNul = TRUE
    )
    names(df) <- c("s", "f")
    ScatFac[[length(ScatFac) + 1]] <- stats::smooth.spline(df$s, df$f)
  }

  names(ScatFac) <- c("H", "He", "Li", "Li+1", "Be", "Be+2", "B", "B+2", "B+3", "C", "N", "O", "O-1", "O-2", "O-1.9", "O-1.4", "F", "F-1", "Ne", "Na", "Na+1", "Mg", "Mg+1", "Mg+2", "Mg+1.9", "Al", "Al+1", "Al+2", "Al+3", "Si", "Si+1", "Si+2", "Si+3", "Si+4", "P", "S", "S-1", "Cl", "Cl-1", "Ar", "K", "K+1", "Ca", "Ca+1", "Ca+2", "Sc", "Sc+1", "Sc+3", "Ti", "Ti+2", "Ti+3", "Ti+4", "V", "V+2", "V+3", "V+5", "Cr", "Cr+2", "Cr+3", "Mn", "Mn+1", "Mn+2", "Mn+3", "Mn+4", "Fe", "Fe+1", "Fe+2", "Fe+3", "Co", "Co+1", "Co+2", "Co+3", "Co+1.2", "Ni", "Ni+1", "Ni+2", "Ni+3", "Ni+0.7", "Cu", "Cu+1", "Cu+2", "Zn", "Zn+2", "Ga", "Ga+3", "Ge", "Ge+4", "As", "Se", "Br", "Br-1", "Kr", "Rb", "Rb+1", "Sr", "Sr+2", "Y", "Y+3", "Zr", "Zr+4", "Nb", "Nb+3", "Nb+5", "Mo", "Mo+3", "Mo+5", "Tc", "Ru", "Rh", "Pd", "Ag", "Ag+1", "Cd", "Cd+2", "In", "In+3", "Sn", "Sn+2", "Sb", "Sb+3", "Te", "I", "I-1", "Xe", "Cs", "Cs+1", "Ba", "Ba+2", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Gd+3", "Tb", "Dy", "Ho", "Er", "Er+3", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "W+6", "Re", "Os", "Ir", "Pt", "Pt+2", "Pt+4", "Au", "Au+1", "Au+3", "Hg", "Hg+1", "Hg+2", "Tl", "Pb", "Pb+2", "Bi", "Bi+3", "U")

  if (save == TRUE) {
    scfile <- paste0(
      "rgl.cry.dp.demo.Sc.",
      format(Sys.time(), "%Y-%m-%d_%H%M%S.Rda")
    )
    save(ScatFac, file = scfile)
  }
}

## from [1, eq. 3.9]
##
##   f(θ) = (1+E0/m0c^2)/(8π^2 a0)  (λ/sin(θ/2))^2 (Z-fx)
##
## where
##
##   E0
##   m0     9.109e-31 kg, from [1, table 1.1]
##   c      2.998e8 m/s, from [1, table 1.1]
##   m0c^2  511 kev, from [1, table 1.1]
##   a0     0.0529 nm,  from [1, eq. 3.5]
##   fx
##
##
##
##  200 keV ->
##
##   (1+E0/m0c^2)/(8π^2 a0) * (λ/sin(θ/2))^2 (Z-fx)
##   (1 + 200e3/511e3)/(8*pi^2 *0.0529e-9) * (0.1541833e-9/sin(θ/2))^2
##
##
## https://physics.nist.gov/PhysRefData/FFast/html/form.html
##

