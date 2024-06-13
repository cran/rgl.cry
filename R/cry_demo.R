# This file is part of the rgl.cry package
#
# Functions and variables:

#' Examples of using the cry and rgl packages together.
#'
#' Read a file in CIF formats, set the parameters, calculates them, and draws
#' the crystal structure with an axis widget.
#'
#' If no file argument is provided, and `dp_demo()` has been opened without
#' paired `cry_demo()`, the CIF parameters of already opened `dp_demo()` will be
#' used.
#'
#' @param file Optional file in CIF formats. The file can also be specified by
#' URL.
#' @param rf A positive value indicating the scale factor of atom radius.
#' @param type A style of atom displaying such like ball, fill and ball-stick
#' but ball-stick is not implemented.
#' @param zoom A positive value indicating the current scene magnification.
#'
#' @return An integer the device number of the currently window.
#'
#' @export
#'
#' @examples
#' cry_demo()
#' cry_demo(system.file("orthorhombic_p.cif", package = "rgl.cry"))
#'
#' \donttest{
#' if (interactive()) {
#'  cry_demo(file, type = "fill", zoom = 0.5)
#'  cry_demo("https://www.crystallography.net/cod/foo.cif")
#' }
#' }
cry_demo <- function(file = NULL, rf = 1, type = "b", zoom = 1) {
  list(file = file, rf = rf, type = type, zoom = zoom)

  ##
  num <- NULL # Number of labeled atom.
  sites <- NULL # Coordinate of labeled atom.
  pos <- NULL # Positions in Cartesian coordinates of labeled atoms.
  mat01 <- NULL # Crystal to Cartesian transform matrix.
  conn <- NULL # Distance between each atom.

  ## The shape of atoms.
  if (type == "fill") {
    radius.factor <- 1
  } else {
    radius.factor <- 0.3
  }


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

  ##
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
      ## Add a column of element name because the label is not an element name.
      ## It seems to work correctly in the following tests.
      ## gsub("[[:digit:][:punct:]]+", "", "O1")
      ## gsub("[[:digit:][:punct:]]+", "", "O1-")
      ## gsub("[[:digit:][:punct:]]+", "", "Cu2+")
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
  ## If there is only one row, the output is a vector, so a workaround is
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

  ## Translating in 26 directions to obtain surrounding coordinates.
  ## dir26 <- unique(t(combn(c(-1, 0, 1, -1, 0, 1, -1, 0, 1), 3)))
  ## dir26 <- dir26[apply(dir26, 1, function(x) any(x != 0)), ] # remove 0, 0, 0
  dir26 <- rbind(
    c(-1, 0, 1), c(-1, 0, -1), c(-1, 0, 0), c(-1, 1, -1),
    c(-1, 1, 0), c(-1, 1, 1), c(-1, -1, 0), c(-1, -1, 1),
    c(-1, -1, -1), c(0, 1, -1), c(0, 1, 0), c(0, 1, 1),
    c(0, -1, 0), c(0, -1, 1), c(0, -1, -1), c(0, 0, 1),
    c(0, 0, -1), c(1, -1, 0), c(1, -1, 1), c(1, -1, -1),
    c(1, 0, 1), c(1, 0, -1), c(1, 0, 0), c(1, 1, -1),
    c(1, 1, 0), c(1, 1, 1)
  )
  for (i in seq_len(num)) {
    tmp <- matrix(apply(dir26, 1, function(v) {
      apply(pos.cry[[i]], 1, function(w) v + w)
    }), ncol = 3, byrow = TRUE)

    tmp <- unique(tmp)
    pos.cry[[i]] <- rbind(pos.cry[[i]], tmp)
  }

  ## Expands the ROI by the specified percentage.
  e <- 0
  for (i in seq_len(num)) {
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[1] >= 0 - e)))
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[1] <= 1 + e)))
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[2] >= 0 - e)))
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[2] <= 1 + e)))
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[3] >= 0 - e)))
    pos.cry[[i]] <-
      subset(pos.cry[[i]], apply(pos.cry[[i]], 1, function(r) all(r[3] <= 1 + e)))
  }

  ## remove duplication
  for (i in seq_len(num)) pos.cry[[i]] <- unique(pos.cry[[i]])


  ## Transform to Cartesian coordinates and store in pos.
  ##
  ## https://en.wikipedia.org/wiki/Fractional_coordinates
  ##
  ##   mat01 <- matrix(c(
  ##     a, b*Cg, c*Cb,
  ##     0, b*Sg, c*(Ca-Cb*Cg)/Sg,
  ##     0, 0, c*sqrt(1-Ca^2-Cb^2-Cg^2+2*Ca*Cb*Cg)/Sg
  ##   ), nrow = 3, byrow = TRUE)
  ##

  ## mat01 <- cry::xtal_mat01(a, b, c, aa, bb, cc)
  mat01 <- cry::xtal_mat02(a, b, c, aa, bb, cc)

  pos <- NULL
  for (i in seq_len(num)) {
    pos[[i]] <- pos.cry[[i]]
    pos[[i]] <- t(apply(pos[[i]], 1, function(v) v %*% mat01))
  }


  ## Calculate distance between each atom.
  conn <- measureDistance(pos)



  ## ------------------------------------------------------------
  ## Define a drawing function.

  drawCry <- function(type) {
    umat <- rgl::par3d("userMatrix")

    ## Atoms
    for (i in seq_len(num)) {
      rgl::spheres3d(pos[[i]],
        r = radius.factor * elemProp[sites[i, 5], ]$radius,
        color = elemProp[sites[i, 5], ]$color,
        tag = "atom"
      )
    }

    ## Crystal coordinates axis and label
    oo <- c(0, 0, 0)
    a1 <- as.vector(t(c(1, 0, 0) %*% mat01))
    a2 <- as.vector(t(c(0, 1, 0) %*% mat01))
    a3 <- as.vector(t(c(0, 0, 1) %*% mat01))

    lines <- rbind(
      oo, a1, oo, a2, (a1 + a2), a1, (a1 + a2), a2,
      a3, (a1 + a3), a3, (a2 + a3), (a1 + a2 + a3), (a1 + a3), (a1 + a2 + a3), (a2 + a3),
      oo, a3, a1, (a1 + a3), a2, (a2 + a3), (a1 + a2), (a1 + a2 + a3)
    )
    rgl::segments3d(lines, col = "grey")



    ## Obtain the combination that has the minimum distance between the atom of
    ## elem[[1]] and the atom of elem[[2]].

    return() # not currently fully runnable.

    smallest <- min(conn[[1]][[2]])
    indices <- which(conn[[1]][[2]] == smallest)
    sticks <- arrayInd(indices, dim(conn[[1]][[2]]))
    sticks <- sticks[, c(2, 1)] # Swapping columns

    ## Draw it as cylinder.
    for (i in seq(1, dim(sticks)[1], by = 1)) {
      c <- rgl::cylinder3d(
        rbind(
          pos[[1]][sticks[i, 1], ],
          pos[[2]][sticks[i, 2], ]
        ),
        radius = 0.05, side = 10
      )
      rgl::shade3d(rgl::addNormals(c), col = "red")
    }
  }




  ## ------------------------------------------------------------

  ## Open device.
  cry.dev <- rgl::open3d()

  inst <- pkg$inst
  if (is.na(inst[nrow(inst), "cry.dev"])) {
    inst[nrow(inst), "cry.dev"] <- cry.dev
  } else {
    inst[nrow(inst) + 1, "cry.dev"] <- cry.dev
  }


  ## Delete child scence and reset the callback to the default.
  try(
    {
      if (length(rgl::subsceneInfo()$parent) != 0) {
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
  cry.root.id <- rgl::currentSubscene3d() # subsceneInfo()$id
  ## 1 rgl::par3d(mouseMode = c("none", "none", "none", "none", "none"))



  ## ------------------------------------------------------------
  ## Subscene for widget

  rgl::newSubscene3d(newviewport = c(0, 0, 150, 150), mouseMode = "replace")

  ## Save the scene id and set the mouseMode for this scene.
  cry.widget.id <- rgl::currentSubscene3d() # subsceneInfo()$id
  ## 1  rgl::par3d(mouseMode = c("none", "none", "none", "none", "none"))

  ## Place the dummy sphere to prevent the draw area from modification.
  ## r=0 is prevention of protrusion.
  rgl::spheres3d(100*frame, r = 0, color = "green", alpha = 0) #
  rgl::par3d(zoom = 0.2) #

  ## Axis widget
  oo <- c(0, 0, 0)
  a1 <- as.vector(t(c(1, 0, 0) %*% mat01))
  a2 <- as.vector(t(c(0, 1, 0) %*% mat01))
  a3 <- as.vector(t(c(0, 0, 1) %*% mat01))

  lines <- rbind(oo, a1, oo, a2, oo, a3)
  rgl::segments3d(lines,
    col = rbind("red", "red", "green", "green", "blue", "blue"),
    lwd = 2
  )

  a0 <- as.vector(t(c(2, 0, 0) %*% mat01))
  b0 <- as.vector(t(c(0, 2, 0) %*% mat01))
  c0 <- as.vector(t(c(0, 0, 2) %*% mat01))
  rgl::text3d(a0, texts = "a", cex = 1.0, col = "black")
  rgl::text3d(b0, texts = "b", cex = 1.0, col = "black")
  rgl::text3d(c0, texts = "c", cex = 1.0, col = "black")


  ## ------------------------------------------------------------
  ## Subscene on the top of scenes for event handling.

  rgl::newSubscene3d(newviewport = c(0, 0, 500, 500), mouseMode = "replace")

  ## Save the scene id and set the mouseMode for this scene.
  cry.panel.id <- rgl::currentSubscene3d()
  ## 1 rgl::par3d(mouseMode = c("none", "user", "none", "none", "pull"))



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
    ##   cry.dev, cry.widget.id, cry.root.id and cry.panel.id

    rgl::set3d(cry.dev, silent = TRUE)
    rgl::useSubscene3d(cry.panel.id)
    umat <- rgl::par3d("userMatrix") # get for


    ## Retrieve the cry and dp pair of this instance.
    inst <- pkg$inst # Get the current list of instance.
    idx <- which(inst$cry.dev == cry.dev)
    start$dp.dev <<- inst[idx, "dp.dev"]
    start$dp.widget.id <<- inst[idx, "dp.widget.id"]
    start$dp.root.id <<- inst[idx, "dp.root.id"]
    start$dp.panel.id <<- inst[idx, "dp.panel.id"]


    ## The rotation is reset to its original value when the mouse is
    ## double-clicked.
    if (time.difference > 100 && time.difference < 200) {
      ## rgl::pop3d(tag = c("rlpoints0","rlpoints1","rlpoints2"))
      umat <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
      rgl::par3d(subscene = cry.widget.id, userMatrix = umat)
      rgl::par3d(subscene = cry.root.id, userMatrix = umat)
      rgl::par3d(subscene = cry.panel.id, userMatrix = umat)
      drawCry()

      ## Update dp_demo drawings if exists.
      if (!is.na(start$dp.dev)) {
        rgl::set3d(start$dp.dev, silent = TRUE)
        rgl::par3d(subscene = start$dp.widget.id, userMatrix = umat)
        rgl::par3d(subscene = start$dp.root.id, userMatrix = umat)
        rgl::par3d(subscene = start$dp.panel.id, userMatrix = umat)
        rgl::useSubscene3d(start$dp.panel.id)
        idx <- which(inst$dp.dev == start$dp.dev)
        inst[[idx, "drawDp"]]()
        rgl::set3d(cry.dev, silent = TRUE)
      }
    }

    ## Save the current mouse cursor position and umat.
    start$x <<- x
    start$y <<- y
    start$umat <<- umat
  }

  ## Rotate and redraw by draging the mouse.
  update <- function(x, y) {

    ## Get the umat at the begining.
    umat <- start$umat # call begin then call update without sequencially
    viewport <- rgl::par3d("viewport", dev = cry.dev)

    w <- viewport[["width"]]
    h <- viewport[["height"]]
    x <- (x - start$x) / w / 1 # 1 is enpirically derived.
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
    rgl::set3d(cry.dev, silent = TRUE)
    rgl::par3d(subscene = cry.widget.id, userMatrix = umat)
    rgl::par3d(subscene = cry.root.id, userMatrix = umat)
    rgl::par3d(subscene = cry.panel.id, userMatrix = umat)
    rgl::useSubscene3d(cry.panel.id)
    ## drawCry() # not necessary

    ## If a pair of dp_demo exists, perform the same operation.
    if (!is.na(start$dp.dev)) {
      rgl::set3d(start$dp.dev, silent = TRUE)
      rgl::par3d(subscene = start$dp.widget.id, userMatrix = umat)
      rgl::par3d(subscene = start$dp.root.id, userMatrix = umat)
      rgl::par3d(subscene = start$dp.panel.id, userMatrix = umat)
      rgl::useSubscene3d(start$dp.panel.id)
      inst <- pkg$inst # Get the current list of instance.
      idx <- which(inst$dp.dev == start$dp.dev)
      inst[[idx, "drawDp"]]()
      rgl::set3d(cry.dev, silent = TRUE)
    }

  }



  ## ------------------------------------------------------------
  ## Start.

  ## Set the callback for the scene cry.panel.id.
  rgl::useSubscene3d(cry.panel.id)
  rgl::rgl.setMouseCallbacks(1, begin, update)
  ## rgl.setWheelCallback(rotate)

  ## Initial view setting.
  rgl::par3d(zoom = zoom)

  ##
  ## cat("Use this panel ID to change the view settings: ", cry.panel.id)

  ##
  drawCry()

  ## Set package global variables
  inst[[nrow(inst), "uc"]] <- uc
  inst[[nrow(inst), "ruc"]] <- ruc
  inst[[nrow(inst), "lCIF"]] <- lCIF
  inst[nrow(inst), "cry.root.id"] <- cry.root.id
  inst[nrow(inst), "cry.widget.id"] <- cry.widget.id
  inst[nrow(inst), "cry.panel.id"] <- cry.panel.id
  assign("inst", inst, pkg)
  ## assign("drawCry", drawCry, pkg)

  ##
  return(as.numeric(cry.dev))
}


## ------------------------------------------------------------
## Element properties
##
## https://en.wikipedia.org/wiki/Atomic_radius
## https://sciencenotes.org/molecule-atom-colors-cpk-colors/

elemPropRaw <- c(
  "H", "#FFFFFF", 25,
  "He", "#D9FFFF", 31, # cal.
  "Li", "#CC80FF", 145,
  "Be", "#C2FF00", 105,
  "B", "#FFB5B5", 85,
  "C", "#909090", 70,
  "N", "#3050F8", 65,
  "O", "#FF0D0D", 60,
  "F", "#90E050", 50,
  "Ne", "#B3E3F5", 38, # cal.
  "Na", "#AB5CF2", 180,
  "Mg", "#8AFF00", 150,
  "Al", "#BFA6A6", 125,
  "Si", "#F0C8A0", 110,
  "P", "#FF8000", 100,
  "S", "#FFFF30", 100,
  "Cl", "#1FF01F", 100,
  "Ar", "#80D1E3", 71, # cal.
  "K", "#8F40D4", 220,
  "Ca", "#3DFF00", 180,
  "Sc", "#E6E6E6", 160,
  "Ti", "#BFC2C7", 140,
  "V", "#A6A6AB", 135,
  "Cr", "#8A99C7", 140,
  "Mn", "#9C7AC7", 140,
  "Fe", "#E06633", 140,
  "Co", "#F090A0", 135,
  "Ni", "#50D050", 135,
  "Cu", "#C88033", 135,
  "Zn", "#7D80B0", 130,
  "Ga", "#C28F8F", 130,
  "Ge", "#668F8F", 125,
  "As", "#BD80E3", 115,
  "Se", "#FFA100", 115,
  "Br", "#A62929", 115,
  "Kr", "#5CB8D1", 88, # cal.
  "Rb", "#702EB0", 235,
  "Sr", "#00FF00", 200,
  "Y", "#94FFFF", 180,
  "Zr", "#94E0E0", 155,
  "Nb", "#73C2C9", 145,
  "Mo", "#54B5B5", 145,
  "Tc", "#3B9E9E", 135,
  "Ru", "#248F8F", 130,
  "Rh", "#0A7D8C", 135,
  "Pd", "#006985", 140,
  "Ag", "#C0C0C0", 160,
  "Cd", "#FFD98F", 155,
  "In", "#A67573", 155,
  "Sn", "#668080", 145,
  "Sb", "#9E63B5", 145,
  "Te", "#D47A00", 140,
  "I", "#940094", 140,
  "Xe", "#429EB0", 108, # cal.
  "Cs", "#57178F", 260,
  "Ba", "#00C900", 215,
  "La", "#70D4FF", 195,
  "Ce", "#FFFFC7", 185,
  "Pr", "#D9FFC7", 185,
  "Nd", "#C7FFC7", 185,
  "Pm", "#A3FFC7", 185,
  "Sm", "#8FFFC7", 185,
  "Eu", "#61FFC7", 185,
  "Gd", "#45FFC7", 180,
  "Tb", "#30FFC7", 175,
  "Dy", "#1FFFC7", 175,
  "Ho", "#00FF9C", 175,
  "Er", "#00E675", 175,
  "Tm", "#00D452", 175,
  "Yb", "#00BF38", 175,
  "Lu", "#00AB24", 175,
  "Hf", "#4DC2FF", 155,
  "Ta", "#4DA6FF", 145,
  "W", "#2194D6", 135,
  "Re", "#267DAB", 135,
  "Os", "#266696", 130,
  "Ir", "#175487", 135,
  "Pt", "#D0D0E0", 135,
  "Au", "#FFD123", 135,
  "Hg", "#B8B8D0", 150,
  "Tl", "#A6544D", 190,
  "Pb", "#575961", 180,
  "Bi", "#9E4FB5", 160,
  "Po", "#AB5C00", 190,
  "At", "#754F45", 127, # cal.
  "Rn", "#428296", 120, # cal.
  "Fr", "#420066", 10, # unknown
  "Ra", "#007D00", 215,
  "Ac", "#70ABFA", 195,
  "Th", "#00BAFF", 180,
  "Pa", "#00A1FF", 180,
  "U", "#008FFF", 175,
  "Np", "#0080FF", 175,
  "Pu", "#006BFF", 175,
  "Am", "#545CF2", 175,
  "Cm", "#785CE3", 10, # unknown
  "Bk", "#8A4FE3", 10, # unknown
  "Cf", "#A136D4", 10, # unknown
  "Es", "#B31FD4", 10, # unknown
  "Fm", "#B31FBA", 10, # unknown
  "Md", "#B30DA6", 10, # unknown
  "No", "#BD0D87", 10, # unknown
  "Lr", "#C70066", 10, # unknown
  "Rf", "#CC0059", 10, # unknown
  "Db", "#D1004F", 10, # unknown
  "Sg", "#D90045", 10, # unknown
  "Bh", "#E00038", 10, # unknown
  "Hs", "#E6002E", 10, # unknown
  "Mt", "#EB0026", 10 # unknown
)

elemProp <- data.frame(
  color = elemPropRaw[c(FALSE, TRUE, FALSE)],
  radius = as.numeric(elemPropRaw[c(FALSE, FALSE, TRUE)]) / 100, # pm to Å
  stringsAsFactors = FALSE
)
rownames(elemProp) <- elemPropRaw[c(TRUE, FALSE, FALSE)]
## elemProp["U",,]$color
## elemProp["U",]$dadius


## ------------------------------------------------------------
## Distance between each atom.

## For example, one possible result is as follows:
##
## [[1]][[2]]
##          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
## [1,] 2.304494 2.304494 2.304494 2.304494 4.412769 4.412769 4.412769 4.412769
## [2,] 5.799515 2.304494 4.412769 4.412769 6.913481 4.412769 4.412769 2.304494
## [3,] 5.799515 4.412769 2.304494 4.412769 4.412769 2.304494 6.913481 4.412769
## [4,] 5.799515 4.412769 4.412769 2.304494 4.412769 2.304494 4.412769 2.304494
## [5,] 4.412769 4.412769 2.304494 2.304494 2.304494 2.304494 5.799515 4.412769
## [6,] 6.913481 4.412769 4.412769 4.412769 5.799515 2.304494 5.799515 2.304494
## [7,] 4.412769 2.304494 2.304494 4.412769 5.799515 4.412769 5.799515 4.412769
## [8,] 4.412769 2.304494 4.412769 2.304494 5.799515 4.412769 2.304494 2.304494
##
## this means the distance between the 1st and 2nd atom combination and
##   column index specifies 1st atom,
##   row index specifies 2nd atom.

## Filter out and collect conbination we need.

measureDistance <- function(p) {
  pos <- p
  conn <- list()
  for (i in seq_along(pos)) {
    conn[[i]] <- sapply(seq_along(pos), function(j) {
      apply(
        pos[[i]], 1, function(w) {
          apply(pos[[j]], 1, function(v) {
            round(stats::dist(rbind(w, v)), digits = 6)
          })
        } # rounding off a 6
      )
    })
  }
  return(conn)
}
