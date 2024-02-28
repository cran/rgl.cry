## This file is part of the rgl.cry package
##
## Functions and variables:

pkg <- NULL

.onLoad <- function(...) {
  pkg <<- new.env()

  ## The selected objects.
  assign("selection", NULL, pkg)

  ## The function to select and identify atoms and reciprocal lattice points.
  assign("select", NULL, pkg)

  ## Call to draw reciprocal lattice points when the orientation changed.
  assign("drawDp", NULL, pkg)

  ## Call to draw crystal structure when the orientation changed.
  ## assign("drawCry", NULL, pkg)


  ## List of cry_demo and dp_demo pairs which contains device id, scene id, lCIF
  ## output and minimam crystal parameters.  I don't know the mechanism yet, but
  ## Bard taught me that the function's environment is kept alive through a
  ## closure, so we can handle it later.
  ## Don't print(inst) so that the drawDp cause error.
  assign(
    "inst",
    data.frame(
      cry.dev = NA, # Device for cry_demo() of this environment.
      cry.root.id = NA, # Scene id.
      cry.widget.id = NA,
      cry.panel.id = NA,
      dp.dev = NA, # Device for dp_demo() at this environment.
      dp.root.id = NA, # Scene id.
      dp.widget.id = NA,
      dp.panel.id = NA,
      drawDp = I(list(1)), # The draw function for reciprocal lattice.
      lCIF = I(list(1)), # The output of CIF.
      uc = I(list(1)), # The unit cell object.
      ruc = I(list(1)), # The reciprocal unit cell object.
      hkl = I(list(1)) # The set of Miller indices up to a given resolution.
    ),
    pkg
  )
}


.onUnload <- function(...) {
  ## Close all opened devices before unloading.

  inst <- pkg$inst

  dev.list <- rgl::rgl.dev.list()

  dev <- inst$dp.dev[which(!is.na(inst$dp.dev))]
  rgl::close3d(dev = intersect(dev.list, dev), silent = TRUE)

  dev <- inst$cry.dev[which(!is.na(inst$cry.dev))]
  rgl::close3d(dev = intersect(dev.list, dev), silent = TRUE)

  message("The rgl.cry package has been unloaded.")
}
