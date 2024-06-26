#' @title Tour in loon
#' @description An interactive tour in loon
#' @name l_tour
#' @inheritParams loon::l_serialaxes
#' @inheritParams tourr::animate
#' @param numOfTours the number of tours
#' @param as.l_tour return a \code{l_tour} object; see details
#' @param interpolation the steps between two serial projections. The larger the value is,
#' the smoother the transitions would be.
#' @param scaling one of 'variable', 'data', 'observation', 'sphere', or 'none' to specify how the data is scaled.
#' See Details
#' @param group only used for layers. As we scroll the bar, the layers are re-calculated.
#' This argument is used to specify which state is used to set groups (i.e. "color", "linewidth", etc).
#' @param slicing whether to show a sliced scatter plot
#' @param slicingDistance the slicing distance that if the distance between
#' points and the projected plane is less than this distance, points will be
#' preserved; else points will be invisible. The default is \code{NULL} and a suggested
#' value will be given. See details
#' @param envir the \code{\link{environment}} to use.
#' @import tourr methods stats loon tcltk loon.ggplot utils
#' @importFrom tibble tibble
#' @details
#' \itemize{
#' \item {tour_path is a tour generator; available tours are \code{\link{grand_tour}},
#' \code{\link{dependence_tour}}, \code{\link{frozen_tour}}, \code{\link{guided_tour}},
#' \code{\link{planned_tour}}, and etc}
#' \item {Argument \code{as.l_tour}
#' \itemize{
#' \item{If set to \code{TRUE}, the function returns an \code{l_tour} (or an \code{l_tour_compound}) object. Essentially, this object is a list with the first element being a \code{loon} (Tcl) widget and the second element a matrix of projection vectors. The advantage of this setup is that the matrix of projection vectors can be easily accessed using the \code{`[`} function (or the \code{l_cget} function). However, a limitation is that it does not constitute a valid \code{loon} (Tcl) widget-calling \code{l_isLoonWidget} would return \code{FALSE}. Consequently, many of loon's functionalities remain inaccessible.}
#' \item{If set to \code{FALSE}, the function returns either a \code{loon} (Tcl) widget (where calling \code{l_isLoonWidget} would return \code{TRUE}) or an \code{l_compound} object. In this case, the matrix of projection vectors is not directly accessible from it. However, the \code{l_getProjection} function can be used to retrieve an estimated matrix of projection vectors.}
#' }
#' }
#' \item {The \code{scaling} state defines how the data is scaled. The axes display 0 at one end and 1 at the other. For the following explanation assume that the data is in a n x p dimensional matrix. The scaling options are then
#' \tabular{ll}{
#' variable \tab per column scaling\cr
#' observation \tab per row scaling\cr
#' data \tab whole matrix scaling\cr
#' sphere \tab transforming variables to principal components}}
#' \item {The default \code{slidingDistance} is suggested by Laa, U., Cook, D., & Valencia, G. (2020).
#' First, find the maximum Euclidean distance of each observation (centralized), say \code{maxD}.
#' Then, compute the "relative volume" that \code{vRel} = (\code{maxD}^(d - 2))/10, where \code{d}
#' is the dimension of this data set. In the end, the suggested \code{slidingDistance}
#' is given by \code{vRel}^(1/(d - 2))
#' }
#' }
#' @return an \code{l_tour} or an \code{l_tour_compound} object that
#' one can query the \code{loon} states and a matrix projection vectors
#' @seealso \code{\link{l_getProjection}}
#' @export
#' @examples
#' if(interactive() && requireNamespace('tourr')) {
#'   # 2D projection
#'   fl <- tourr::flea[, 1:6]
#'   # different scaling will give very different projections
#'   # in this dataset, scaling 'variable' will give the best separation
#'   p <- l_tour(fl, scaling = 'variable',
#'               color = tourr::flea$species)
#'   l0 <- l_layer_hull(p, group = p["color"],
#'                      color = "red", linewidth = 4)
#'   l1 <- l_layer_density2d(p)
#'   # a `l_tour` object
#'   class(p)
#'
#'   # query the matrix of projection vectors
#'   proj <- p['projection'] # or `l_getProjection(p)`
#'   # suppose the scaling is still 'observation'
#'   new_xy <- as.matrix(
#'     loon::l_getScaledData(data = fl,
#'                           scaling = 'observation')) %*%
#'     proj
#'   plot(new_xy, xlab = "V1", ylab = "V2",
#'        col = loon::hex12tohex6(p['color']))
#'
#'   # A higher dimension projection
#'   # turn the `tour` to 4 dimensional space
#'   s <- l_tour(fl, color = tourr::flea$species,
#'               scaling = "observation",
#'               tour_path = tourr::grand_tour(4L))
#'
#'   # set `as.l_tour` FALSE
#'   p <- l_tour(fl, scaling = 'observation',
#'               color = tourr::flea$species)
#'   class(p)
#'   ## ERROR
#'   ## p["projection"]
#'
#'   # query the estimated matrix of projection vectors
#'   l_getProjection(p)
#'
#'   ##### facet by region
#'   olive <- tourr::olive
#'   p <- with(olive, l_tour(olive[, -c(1, 2)],
#'                           by = region,
#'                           color = area))
#' }
l_tour <- function(data, scaling = c('data', 'variable', 'observation', 'sphere'),
                   by = NULL, on, as.l_tour = TRUE, color = loon::l_getOption("color"),
                   tour_path = tourr::grand_tour(), group = "color", start = NULL,
                   slicing = FALSE, slicingDistance = NULL,
                   numOfTours = 30L, interpolation = 40L, parent = NULL,
                   envir = parent.frame(), ...) {


  stopifnot(
    exprs = {
      is.numeric(numOfTours)
      is.finite(numOfTours)
      is.numeric(interpolation)
      is.finite(interpolation)
    }
  )

  ###### THIS IS A HACK IN `TOURR`! Hopefully, it can be fixed later.
  record <- tibble::tibble(basis = list(), index_val = numeric(),
                           info = character(), method = character(),
                           alpha = numeric(),
                           tries = numeric(), loop = numeric())
  ############################

  scaling <- match.arg(scaling)
  # the dimension of data set
  # data is n * d
  originalData <- data
  # start should be d * k
  originalStart <- start
  n <- nrow(data)

  data <- get_scaledData(originalData, scaling = scaling)
  d <- ncol(data)
  start <- creat_start(data, originalStart, tour_path, d)
  k <- ncol(start)

  env <- environment()

  # projections are
  # l * d * k
  projections <- interpolate_list(data, start = start,
                                  tour_path = tour_path, numOfTours = numOfTours,
                                  interpolation = interpolation)
  # tours are
  # l * ((n * d) * (d * k)) --> l * n * k
  tours <- tour_list(data, projections)

  len_color <- length(color)
  if (len_color > 1) {
    if (len_color != n) {
      color <- rep_len(color, n)
    }
  } else {
    if(is.na(color)) color <- loon::l_getOption("color")
    color <- rep(color, n)
  }

  new.toplevel <- FALSE
  if(is.null(parent)) {
    new.toplevel <- TRUE
    # create parent
    parent <- loon::l_toplevel()
  }

  if(is.null(by)) {

    byDataFrame <- NULL

    # in single widget
    # parent is converted to child
    subwin <- loon::l_subwin(parent, "tour")
    child <- as.character(tcltk::tcl('frame', subwin))
    tcltk::tktitle(parent) <- paste("Tour", "--path:", child)
    parent <- child
  } else {

    # in compound widget
    # parent is the tk top level
    # child is the valid path
    dotArgs <- list(...)
    dotArgs$color <- color

    l_className <- if(d == 1) {
      "l_hist"
    } else if(d == 2) {
      "l_plot"
    } else {
      "l_serialaxes"
    }
    byDataFrame <- by2Data(by, on, bySubstitute = substitute(by), n = n,
                           args = dotArgs, l_className = l_className)

    tours <- lapply(tours,
                    function(tour) {
                      split(as.data.frame(tour),
                            f = as.list(byDataFrame),
                            drop = FALSE)
                    })

    if(!inherits(by, "formula"))
      by <- byDataFrame
    child <- paste0(loon::l_subwin(parent, ""), "facet")
    tcltk::tktitle(parent) <- paste("Tour", "--path:", child)
  }

  # p could be a `l_plot`, `l_hist`, `l_serialaxes`
  # in addition, it could be a `l_facet` object as well.
  p <- l_creat_tour_widget(data, scaling = "none",
                           color = color, tour_path = tour_path,
                           start = start,  parent = parent, by = by, on, ...)

  axes <- l_layer_getAxes(p)
  labels <- l_layer_getLabels(p)
  ################################ tour tools ################################
  # Set up position
  from <- 0
  var <- 0
  varTcl <- tcltk::tclVar(var)
  # `length(projections)` is less and equal than numOfTours * interpolation
  to <- length(projections)
  # tour slider bar
  tourBar <- as.character(
    tcltk::tcl('scale',
               as.character(loon::l_subwin(child,'scale')),
               orient = 'vertical', variable = varTcl,
               showvalue = 0, from = from, to = to,
               resolution = 1)
  )

  # refresh button
  refreshButton <- as.character(
    tcltk::tcl('button',
               as.character(loon::l_subwin(child,'refresh button')),
               text = "refresh",
               bg = "grey80",
               fg = "black",
               borderwidth = 2,
               relief = "raised"))

  # scale radio button
  scalingVar <- tcltk::tclVar(scaling)
  scaleRadioButtons <- sapply(
    as.character(formals(l_tour)[["scaling"]])[-1],
    function(scale) {
      as.character(
        tcltk::tcl('radiobutton',
                   as.character(loon::l_subwin(child,
                                               paste0(c('radiobutton', scale),
                                                      collapse = "-"))
                   ),
                   text = scale,
                   variable = scalingVar,
                   value = scale))
    })

  onStepUpButton <- tryCatch(
    expr = {
      path <- file.path(find.package(package = 'loon.tourr'), "images")
      upbutton <- tcltk::tkimage.create("photo",  tcltk::tclVar(),
                                        file=paste0(path, "/up.png"))

      # one step up
      as.character(
        tcltk::tcl('button',
                   as.character(loon::l_subwin(child, 'one step up button')),
                   image = upbutton,
                   bg = "grey92",
                   borderwidth = 1,
                   relief = "raised"))

    },
    error = function(e) {
      # one step up
      as.character(
        tcltk::tcl('button',
                   as.character(loon::l_subwin(child, 'one step up button')),
                   text = "up",
                   bg = "grey92",
                   borderwidth = 1,
                   relief = "raised"))
    }
  )

  # one step down
  onStepDownButton <- tryCatch(
    expr = {
      path <- file.path(find.package(package = 'loon.tourr'), "images")
      downbutton <- tcltk::tkimage.create("photo",  tcltk::tclVar(),
                                          file=paste0(path, "/down.png"))
      as.character(
        tcltk::tcl('button',
                   as.character(loon::l_subwin(child, 'one step down button')),
                   image = downbutton,
                   bg = "grey92",
                   borderwidth = 1,
                   relief = "raised"))
    },
    error = function(e) {
      as.character(
        tcltk::tcl('button',
                   as.character(loon::l_subwin(child, 'one step down button')),
                   text = "down",
                   bg = "grey92",
                   borderwidth = 1,
                   relief = "raised"))
    }
  )

  rowcols <- tk_get_row_and_columns(p, by = by, title = list(...)$title)

  tk_grid_pack_tools(tourBar = tourBar, refreshButton = refreshButton,
                     scaleRadioButtons = scaleRadioButtons,
                     onStepUpButton = onStepUpButton,
                     onStepDownButton = onStepDownButton,
                     parent = child,
                     pack = new.toplevel,
                     row = rowcols$row,
                     column = rowcols$column,
                     by = by)

  # the setting of the color is used in the layer (`l_layer`) objects
  # if multiple color is displayed
  # as we pack layers, the layer will be grouped by color
  if(!is.null(byDataFrame))
    color <- split(color, f = as.list(byDataFrame), drop = FALSE)
  # initial settings
  axesLength <- 0.2
  count <- 0
  countOld <- 0
  # `varOld` is only used in trail layer
  varOld <- 0
  scalingOld <- scaling
  # step <- 0
  toOld <- to

  # group variable
  group <- tryCatch(
    expr = {
      p[group]
    },
    error = function(e) ""
  )
  initialTour <- l_initialTour(p)

  tour <- if(as.l_tour) {
    if(is(p, "l_compound")) {
      structure(
        stats::setNames(list(p, start),
                        c(as.character(child), "projection")),
        class = c("l_tour_compound", "loon"))
    } else {
      structure(
        stats::setNames(list(p, start),
                        c(as.character(child), "projection")),
        class = c("l_tour", "loon"))
    }
  } else {
    p
  }

  states <- l_getNDimStates(p)
  statesNames <- names(states)

  if(slicing) {

    slicingDistance <- get_slicingDistance(slicingDistance,
                                           originalData, d)
    suppressMessages(
      loon::l_bind_state(p, event = c("glyph", "color", "selected",
                                      "active", "size"),
                         callback = function(e) {

                           original.linkingKey <- states$linkingKey
                           current.linkingKey <- p['linkingKey']

                           states <<- stats::setNames(
                             lapply(statesNames,
                                    function(name) {

                                      s <- states[[name]]
                                      if(name == "data") {
                                        s[original.linkingKey %in% current.linkingKey, ] <- p[name]
                                      } else {
                                        s[original.linkingKey %in% current.linkingKey] <- p[name]
                                      }
                                      s
                                    }),
                             statesNames)
                         })
    )
  }

  l_object_name <- character(0L)

  update <- function(...) {

    # callback functions have side effects
    scalingVar <- tcltk::tclvalue(scalingVar)
    # data, start, initialTour, projections
    # tours, scalingOld are modified
    callback_scaling(scalingVar, scalingOld, tour_path,
                     originalData, originalStart, d, k, start,
                     numOfTours, interpolation, by = by, env = env)

    # tours, projections and countOld are modified
    callback_refresh(count, countOld, data, start,
                     tour_path, numOfTours, interpolation,
                     by = by, env = env)

    # as the slider bar is moved,
    # the `varTcl` is updated automatically
    # then, update the `var` manually
    var <<- as.numeric(tcltk::tclvalue(varTcl))

    if(as.l_tour) {
      # A compromised way to pack a projection object in loon
      # Create a list "list(loon, projection)", then modify the class of the list (force it to be a loon widget)
      # Advantage:
      #   1. No need to touch the tcl or change the main structure of loon
      #   2. The projection matrix can be updated as the interface changes.
      # Drawbacks:
      #   1. Have to create a bunch of functions, like `[.l_tour`, `[<-.l_tour`, ... to cover
      #   2. A `l_tour` must be assigned to a real object.
      #      ########### Example
      #      # in your console, if you call
      #      l_tour(iris)
      #      # Instead of
      #      p <- l_tour(iris)
      #      `projection` matrix would not be updated
      #   3. Several objects may point to one single loon widget.

      # update l_object_name
      # find the loon object name in `envir`
      ls_objects <- ls(envir = envir)
      obj_name <- lapply(ls_objects,
                         function(obj) {
                           tourobj <- get(obj, envir = envir)
                           pname <- l_path_name(tourobj, child)
                           if(length(pname) == 0) return(NULL)
                           if(grepl(as.character(child), pname) && (l_isTour(tourobj) || loon::l_isLoonWidget(tourobj) || l_isCompound(tourobj))) return(obj)

                           else return(NULL)
                         })
      assign("l_object_name", unlist(obj_name), envir = env)

      if(var > 0 && !is.null(l_object_name)) {
        tour <- l_configure(
          tour,
          projection = projections[[var]]
        )
        # This loop is very dangerous
        # Reason: more than one objects would point to one single loon widget
        # So far, I do not have any better idea
        # FIXME
        # Find Better ways to extract the matrix of projection vectors (side effect is not suggested)
        lapply(l_object_name,
               function(x) {
                 tourobj <- get(x, envir = envir)
                 pname <- l_path_name(tourobj, child)
                 if(grepl(as.character(child), pname) && (l_isTour(tourobj) || loon::l_isLoonWidget(tourobj) || l_isCompound(tourobj)))
                   assign(x, tour, envir = envir)
                 else {
                   warning("The object `",
                           x,
                           "` is assigned to other object ",
                           "so that the `projections` will not be updated",
                           call. = FALSE,
                           immediate. = TRUE)}
               })
      }
    }

    callback_plot(widget = p, initialTour = initialTour, tours = tours, var = var,
                  data = originalData, start = start, color = color, group = group,
                  varOld = varOld, projections = projections,
                  axesLength = axesLength, axes = axes, labels = labels,
                  states = states, slicing = slicing,
                  slicingDistance = slicingDistance)

    # step <<- step + 1
    varOld <<- var
    to <<- length(projections)

    if(to != toOld) {
      toOld <<- to
      tcltk::tkconfigure(tourBar, to = to)
    }
  }
  # as tourBar or refresh button is modified, update function will be executed
  tcltk::tkconfigure(tourBar, command = update)
  tcltk::tkconfigure(refreshButton,
                     command = function(...) {
                       count <<- count + 1
                       update(...)
                     })

  tcltk::tkconfigure(onStepUpButton,
                     command = function(...) {

                       if(var == 0) return(NULL)

                       varOld <<- var
                       var <<- var - 1
                       varTcl <<- tcltk::tclVar(var)
                       tcltk::tkconfigure(tourBar,
                                          variable = varTcl)
                       update(...)
                     })


  tcltk::tkconfigure(onStepDownButton,
                     command = function(...) {

                       if(var == to) return(NULL)

                       varOld <<- var
                       var <<- var + 1
                       varTcl <<- tcltk::tclVar(var)
                       tcltk::tkconfigure(tourBar,
                                          variable = varTcl)
                       update(...)
                     })

  lapply(scaleRadioButtons,
         function(scale_radio_button) {
           tcltk::tkconfigure(scale_radio_button,
                              command = update)
         })

  return(tour)
}
