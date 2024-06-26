% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l_layer_density2d.R
\name{l_layer_density2d}
\alias{l_layer_density2d}
\title{2D density layer}
\usage{
l_layer_density2d(
  widget,
  x,
  y,
  h,
  n = 25L,
  lims = NULL,
  color = "black",
  linewidth = 1,
  nlevels = 10,
  levels = NULL,
  label = "density2d",
  parent = "root",
  index = 0,
  group = NULL,
  active = TRUE,
  ...
)
}
\arguments{
\item{widget}{`loon` widget path name as a string}

\item{x}{The coordinates of \code{x}. See details}

\item{y}{The coordinates of \code{y}. See details}

\item{h}{
vector of bandwidths for x and y directions.  Defaults to
normal reference bandwidth (see \code{\link[MASS]{bandwidth.nrd}}). A scalar
value will be taken to apply to both directions.
}

\item{n}{
Number of grid points in each direction.  Can be scalar or a length-2
integer vector.
}

\item{lims}{
The limits of the rectangle covered by the grid as \code{c(xl, xu, yl, yu)}.
}

\item{color}{color of each contour}

\item{linewidth}{the line width}

\item{nlevels}{As described in \code{grDevices::contourLines}: number of
contour levels desired iff \code{levels} is not supplied.}

\item{levels}{As described in \code{grDevices::contourLines}: numeric vector
of levels at which to draw contour lines.}

\item{label}{label used in the layers inspector}

\item{parent}{parent group layer}

\item{index}{of the newly added layer in its parent group}

\item{group}{separate x vector or y vector into a list by group}

\item{active}{a logical determining whether points appear or not (default is \code{TRUE} for all points).
If a logical vector is given of length equal to the number of points,
then it identifies which points appear (\code{TRUE}) and which do not (\code{FALSE}).}

\item{...}{other arguments to modify \code{l_layer_line}.}
}
\value{
an \code{l_layer} widget
}
\description{
Two-dimensional kernel density estimation with an axis-aligned bivariate normal kernel
}
\examples{
if(interactive()) {
p <- l_plot(iris, color = iris$Species)
l <- l_layer_density2d(p)
}
}
