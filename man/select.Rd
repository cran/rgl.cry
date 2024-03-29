% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select}
\alias{select}
\title{Select atoms or reciprocal lattice points.}
\usage{
select(dev = NULL, verbose = TRUE)
}
\arguments{
\item{dev}{RGL device to apply.  Defaults to current device.}

\item{verbose}{logical: Should the report be suppressed?}
}
\value{
List of Miller indices or element labels.
}
\description{
Select one or more atoms or reciprocal lattice points in the window.  The
labels and Miller indices of the selected atoms or lattice points will be
displayed.
}
\details{
Selecting atoms or lattice points in the window will include all
z-coordinates.  If you do not want to include all z-coordinates, you will
need to modify the code.
}
\examples{
\donttest{
if (interactive()) {
 select()
 select(dev = 123)
}
}
}
