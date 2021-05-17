## ----setup, include=FALSE, warning=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE,
                      message = FALSE,
                      fig.align = "center", 
                      fig.width = 6, 
                      fig.height = 5,
                      out.width = "80%", 
                      collapse = TRUE,
                      comment = "#>",
                      tidy.opts = list(width.cutoff = 65),
                      tidy = FALSE)
library(knitr)
library(magrittr)
library(loon.tourr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(MASS)
set.seed(12314159)
imageDirectory <- "./images/tours"
dataDirectory <- "./data/tours"
path_concat <- function(path1, ..., sep="/") {
  # The "/" is standard unix directory separator and so will
  # work on Macs and Linux.
  # In windows the separator might have to be sep = "\" or 
  # even sep = "\\" or possibly something else. 
  paste(path1, ..., sep = sep)
}

## ----library loon.tourr, warning=FALSE, message=FALSE, error=FALSE------------
library(loon.tourr)

## ----data, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
library(MASS, quietly = TRUE)
kable(head(crabs, 6))

## ----crabs 2D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  color <- rep("skyblue", nrow(crabs))
#  color[crabs$sp == "O"] <- "orange"
#  cr <- crabs[, c("FL", "RW", "CL", "CW", "BD")]
#  p0 <- l_tour(cr, color = color)

## ----crabs 2D gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 1: 2D grand tour"----
include_graphics(path_concat(imageDirectory, "crab2D.gif"))

## ----crabs 4D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  p1 <- l_tour(cr,
#               tour_path = grand_tour(4L),
#               color = color)

## ----crabs 4D gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 2: 4D grand tour"----
include_graphics(path_concat(imageDirectory, "crab4D.gif"))

## ----`l_tour` class, eval = FALSE-----------------------------------------------
#  class(p1) # class(p0)
#  # > [1] "l_tour" "loon"

## ----`l_tour` getLoon, eval = FALSE---------------------------------------------
#  w <- l_getPlots(p1)
#  class(w)
#  # > [1] "l_serialaxes" "loon"

## ----query matrix of projection vectors, eval = FALSE-------------------------
#  round(p1["projection"], 2)
#  # >
#  #       [,1]  [,2]  [,3]  [,4]
#  # [1,] -0.12 -0.93 -0.01 -0.35
#  # [2,] -0.57  0.21  0.63 -0.36
#  # [3,] -0.06  0.14  0.14 -0.38
#  # [4,] -0.37  0.23 -0.76 -0.47
#  # [5,] -0.72 -0.14 -0.12  0.62

## ----parallel, eval = FALSE---------------------------------------------------
#  p1["axesLayout"] <- "parallel"

## ----andrews, eval = FALSE----------------------------------------------------
#  p1["andrews"] <- TRUE

## ----static grid, eval = FALSE------------------------------------------------
#  plot(p1)

## ----static ggplot, eval = FALSE----------------------------------------------
#  loon.ggplot::loon.ggplot(p1)

## ----crabs andrews, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "Figure 1: Andrews curve"----
include_graphics(path_concat(imageDirectory, "andrews.png"))

## ----gt, eval = FALSE---------------------------------------------------------
#  # Default, 2D grand tour
#  pg <- l_tour(cr, tour_path = grand_tour(2L))

## ----holes, eval = FALSE------------------------------------------------------
#  # 2D holes projection pursuit indexes
#  pp_holes <- l_tour(cr, tour_path = guided_tour(holes(), 2L))

## ----cmass, eval = FALSE------------------------------------------------------
#  # 2D CM projection pursuit indexes
#  pp_CM <- l_tour(cr, tour_path = guided_tour(cmass(), 2L))

## ----PCA, eval = FALSE--------------------------------------------------------
#  # 2D LDA projection pursuit indexes
#  pp_LDA <- l_tour(cr,
#                   color = crabs$sex,
#                   tour_path = guided_tour(lda_pp(crabs$sex), 2L))

## ----crabs facets, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  pf <- l_tour(cr,
#               by = crabs$sex,
#               color = color)

## ----crabs facets gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 3: facets"----
include_graphics(path_concat(imageDirectory, "facets.gif"))

## ----crabs pairs, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  pp <- l_tour_pairs(cr,
#                     tour_path = grand_tour(4L),
#                     color = color,
#                     showSerialAxes = TRUE)

## ----crabs pairs gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 4: pairs plot"----
include_graphics(path_concat(imageDirectory, "pairs.gif"))

## ----`l_tour_compound` class, eval = FALSE--------------------------------------
#  class(pp) # or class(pf)
#  # > [1] "l_tour_compound" "loon"

## ----pairs objects, eval = FALSE----------------------------------------------
#  wp <- l_getPlots(pp)
#  wp
#  # >
#  # $x2y1
#  # [1] ".l0.pairs.plot"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $x3y1
#  # [1] ".l0.pairs.plot1"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $x4y1
#  # [1] ".l0.pairs.plot2"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $x3y2
#  # [1] ".l0.pairs.plot3"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $x4y2
#  # [1] ".l0.pairs.plot4"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $x4y3
#  # [1] ".l0.pairs.plot5"
#  # attr(,"class")
#  # [1] "l_plot" "loon"
#  #
#  # $serialAxes
#  # [1] ".l0.pairs.serialaxes"
#  # attr(,"class")
#  # [1] "l_serialaxes" "loon"
#  #
#  # attr(,"class")
#  # [1] "l_pairs"    "l_compound" "loon"

## ----crabs hull, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  # pack layer on top of `p0`
#  l0 <- l_layer_hull(p0, group = crabs$sp)

## ----hull gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 5: layer hull"----
include_graphics(path_concat(imageDirectory, "hull.gif"))

## ----crabs density2D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  # hide the hull
#  l_layer_hide(l0)
#  # density2D
#  l1 <- l_layer_density2d(p0)

## ----density 2D gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 5: layer density"----
include_graphics(path_concat(imageDirectory, "density2D.gif"))

## ----crabs trails, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  # hide the density2D
#  l_layer_hide(l1)
#  # density2D
#  l2 <- l_layer_trails(p0)

## ----crabs trails gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 5: layer density"----
include_graphics(path_concat(imageDirectory, "trails.gif"))

## ----non-interactive layer, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  allx <- unlist(pf['x'])
#  ally <- unlist(pf['y'])
#  layers <- lapply(l_getPlots(pf),
#                   function(p) {
#                     l <- loon::l_layer_points(p,
#                                               x = allx,
#                                               y = ally,
#                                               color = "grey80",
#                                               label = "background")
#                     # set the layer as the background
#                     loon::l_layer_lower(p, l)
#                   })

## ----non-interactive layer gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "Figure 2: non-interactive layer"----
include_graphics(path_concat(imageDirectory, "tour_layer_non_interactive.PNG"))

## ----interactive layer, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  l_layer_callback.background <- function(target, layer, ...) {
#  
#    widget <- l_getPlots(target)
#    layer <- loon::l_create_handle(c(widget, layer))
#  
#    args <- list(...)
#    # the overall tour paths
#    allTours <- args$allTours
#    # the scale bar variable
#    var <- args$var
#    # the current projection (bind both facets)
#    proj <- do.call(rbind, allTours[[var]])
#  
#    loon::l_configure(layer,
#                      x = proj[, 1],
#                      y = proj[, 2])
#  }

## ----interactive layer gif, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", fig.cap = "GIF 7: interactive layer"----
include_graphics(path_concat(imageDirectory, "tour_layer_interactive.gif"))

