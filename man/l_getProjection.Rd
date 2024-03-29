% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l_getProjection.R
\name{l_getProjection}
\alias{l_getProjection}
\title{Query the matrix of projection vectors}
\usage{
l_getProjection(target, data)
}
\arguments{
\item{target}{A object returned by \code{l_tour}}

\item{data}{Original data set}
}
\value{
a matrix of projection vectors
}
\description{
Query the matrix of projection vectors
}
\examples{
if(interactive()) {
  dat <- iris[,-5]
  p <- l_tour(dat, color = iris$Species,
              as.l_tour = FALSE)
  # scroll the bar
  proj <- l_getProjection(p, dat)
  projected_object <- as.matrix(dat) \%*\% proj
  # it will not be exactly the same
  plot(projected_object[,1], projected_object[,2],
       col = hex12tohex6(p['color']))
}
}
