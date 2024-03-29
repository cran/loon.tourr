% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l_tour_pairs.R
\name{tour_pairs}
\alias{tour_pairs}
\alias{l_tour_pairs}
\title{Tour Pairs Plot}
\usage{
l_tour_pairs(
  data,
  scaling = c("data", "variable", "observation", "sphere"),
  tour_path = tourr::grand_tour(4L),
  numOfTours = 30L,
  interpolation = 40L,
  as.l_tour = TRUE,
  connectedScales = c("none", "cross"),
  linkingGroup,
  linkingKey,
  showItemLabels = TRUE,
  itemLabel,
  showHistograms = FALSE,
  histLocation = c("edge", "diag"),
  histHeightProp = 1,
  histArgs = list(),
  showSerialAxes = FALSE,
  serialAxesArgs = list(),
  color = "grey60",
  group = "color",
  start = NULL,
  parent = NULL,
  span = 10L,
  envir = parent.frame(),
  ...
)
}
\arguments{
\item{data}{a data frame with numerical data only}

\item{scaling}{one of 'variable', 'data', 'observation', 'sphere', or 'none' to specify how the data is scaled.
See Details}

\item{tour_path}{tour path generator, defaults to 2d grand tour}

\item{numOfTours}{the number of tours}

\item{interpolation}{the steps between two serial projections. The larger the value is,
the smoother the transitions would be.}

\item{as.l_tour}{return a \code{l_tour} object; see details}

\item{connectedScales}{Determines how the scales of the panels are to be connected.
\itemize{
\item{\code{"cross"}: only the scales in the same row and the same column are connected;}
\item{\code{"none"}: neither "x" nor "y" scales are connected in any panels.}
}}

\item{linkingGroup}{string giving the linkingGroup for all plots. If missing,
a default \code{linkingGroup} will be determined from deparsing the \code{data}.}

\item{linkingKey}{a vector of strings to provide a linking identity for each row of the
\code{data} data.frame.  If missing, a default \code{linkingKey} will be \code{0:(nrows(data)-1)}.}

\item{showItemLabels}{TRUE, logical indicating whether its itemLabel pops up over a point when
the mouse hovers over it.}

\item{itemLabel}{a vector of strings to be used as pop up information when the mouse hovers
over a point.  If missing, the default \code{itemLabel} will be the \code{row.names(data)}.}

\item{showHistograms}{logical (default FALSE) to show histograms of each variable
or not}

\item{histLocation}{one "edge" or "diag", when showHistograms = TRUE}

\item{histHeightProp}{a positive number giving the height of the histograms as a
proportion of the height of the scatterplots}

\item{histArgs}{additional arguments to modify the `l_hist` states}

\item{showSerialAxes}{logical (default FALSE) indication of whether to show a serial axes plot
in the bottom left of the pairs plot (or not)}

\item{serialAxesArgs}{additional arguments to modify the `l_serialaxes` states}

\item{color}{vector with line colors. Default is given by \code{\link[loon]{l_getOption}("color")}.}

\item{group}{only used for layers. As we scroll the bar, the layers are re-calculated.
This argument is used to specify which state is used to set groups (i.e. "color", "linewidth", etc).}

\item{start}{projection to start at, if not specified, uses default
associated with tour path}

\item{parent}{a valid Tk parent widget path. When the parent widget is
specified (i.e. not \code{NULL}) then the plot widget needs to be placed using
some geometry manager like \code{\link{tkpack}} or \code{\link{tkplace}} in
order to be displayed. See the examples below.}

\item{span}{How many column/row occupies for each widget}

\item{envir}{the \code{\link{environment}} to use.}

\item{...}{named arguments to modify the serialaxes states or layouts, see details.}
}
\value{
an \code{l_tour_compound} object that one can query the \code{loon} states
and a matrix projection vectors
}
\description{
A nD tour path with a scatterplot matrix (the default tour is a 4D tour;
by setting `tour_path` to modify the dimension)
}
\examples{
if(interactive() && requireNamespace('tourr')) {
  # q is a `l_pairs` object
  q <- l_tour_pairs(olive[, -c(1:2)],
                    color = olive$region)
  # query the matrix of projection vectors
  proj <- q["projection"]

  # query the `l_compound` widget
  lc <- l_getPlots(q)
  # pack the `density2d` layers
  layer_pack <- lapply(lc, function(w) l_layer_density2d(w))

  #### set `as.l_tour = FALSE`
  # q is a `l_pairs` object
  q <- l_tour_pairs(tourr::flea[, 1:6],
                    as.l_tour = FALSE,
                    color = tourr::flea$species,
                    showHistogram = TRUE,
                    showSerialAxes = TRUE)

  # proj <- q["projection"] # Return a list of `NA`
  # query estimated matrix of projection vectors
  proj <- l_getProjection(q, tourr::flea[, 1:6])
}
}
\seealso{
\code{\link{l_pairs}}, \code{\link{l_tour}}
}
