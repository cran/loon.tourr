## ----setup, include=FALSE, warning=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE,
                      message = FALSE,
                      fig.align = "center", 
                      fig.width = 6, 
                      fig.height = 5,
                      out.width = "60%", 
                      collapse = TRUE,
                      comment = "#>",
                      tidy.opts = list(width.cutoff = 65),
                      tidy = FALSE)
library(knitr)
library(magrittr)
library(loon.tourr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(class, quietly = TRUE)
set.seed(12314159)
imageDirectory <- "./images/classification"
dataDirectory <- "./data/classification"
path_concat <- function(path1, ..., sep="/") {
  # The "/" is standard unix directory separator and so will
  # work on Macs and Linux.
  # In windows the separator might have to be sep = "\" or 
  # even sep = "\\" or possibly something else. 
  paste(path1, ..., sep = sep)
}
library(RDRToolbox)

## ----table, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
data.frame(
  Region = c("North", "South", "Sardinia"),
  Area = c("North-Apulia, South-Apulia, Calabria, Sicily",
           "East-Liguria, West-Liguria, Umbria", 
           "Coastal-Sardinia, Inland-Sardinia")
) %>% 
  kable()

## ----split data, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
set.seed(123)
N <- nrow(olive)
trainId <- sample(seq(N), 
                  size = floor(0.8 * N))
testId <- setdiff(seq(N), trainId)
acids <- setdiff(colnames(olive), c("region", "area"))
trainX <- olive[trainId, acids]
testX <- olive[testId, acids]
trainY <- olive[trainId, "region"]
testY <- olive[testId, "region"]

## ----scaling, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
row.names(trainX) <- NULL
kable(head(trainX))
scalingTrainX <- loon::l_getScaledData(trainX, scaling = "variable")
scalingTestX <- loon::l_getScaledData(testX, scaling = "variable")
kable(head(scalingTrainX), digits = 2)

## ----xgboost, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
knn_pred <- function(trainX, trainY, testX, testY, k = c(5, 10, 20)) {
  len_test <- length(testY)
  vapply(k,
         function(num) {
           yhat <- class::knn(trainX, testX, trainY, k = num)
           sum(yhat == testY)/len_test
         }, numeric(1L))
  
}

low_dim_knn_pred <- function(dims = 2:5, fun, 
                             trainX, trainY, testX, testY, 
                             k = c(5, 10, 20), setNames = TRUE) {
  
  tab <- lapply(dims, 
                function(d) {
                  knn_pred(
                    fun(trainX, d), trainY,
                    fun(testX, d), testY
                  )
                }) %>% 
    as.data.frame() %>%
    as_tibble()
  
  if(setNames) {
    tab <- tab %>%
      setNames(nm = paste0("d = ", dims))
  }
  
  rownames(tab) <- paste0("k = ", k)
  tab
}

## ----p choose k, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
# the number of k
dims <- 2:5
var_names <- colnames(scalingTrainX)
low_dim_names <- c()
K <- ncol(trainX)
pChooseD <- lapply(dims, 
                   function(d) {
                     com <- combn(K, d)
                     pred <- apply(com, 2, 
                                   function(pair) {
                                     knn_pred(trainX[, pair], trainY, 
                                              testX[, pair], testY)
                                   })
                     mean_pred <- apply(pred, 2, mean)
                     id <- which.max(mean_pred)
                     low_dim_names <<- c(low_dim_names, 
                                         paste(var_names[com[, id]], 
                                               collapse = ":"))
                     pred[, id]
                   }) %>% 
  as.data.frame() %>%
  as_tibble() %>%
  setNames(nm = paste0("d = ", dims))
rownames(pChooseD) <- paste0("k = ", c(5, 10, 20))

## ----pairs--------------------------------------------------------------------
names(low_dim_names) <- paste0("d = ", dims)
low_dim_names

## ----nav graph prediction table-----------------------------------------------
kable(pChooseD, row.names = TRUE,
      digits = 3)

## ----PCA, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
trainXPCA <- princomp(scalingTrainX)
testXPCA <- princomp(scalingTestX)
round(trainXPCA$sdev, 2)

## ----PCA knn, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
PCA <- low_dim_knn_pred(2:5, 
                        fun = function(princomp, d) 
                          {princomp$scores[, seq(d)]},
                        trainXPCA,
                        trainY,
                        testXPCA,
                        testY)
kable(PCA, row.names = TRUE,
      digits = 3)

## ----LLE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", eval = FALSE----
#  lle <- low_dim_knn_pred(2:5,
#                          fun = function(data, d) {
#                            LLE(data, dim = d, k = 5)
#                          },
#                          scalingTrainX, trainY,
#                          scalingTestX, testY)
#  kable(lle, row.names = TRUE,
#        digits = 3)

## ----LLE readRDS, echo = FALSE------------------------------------------------
lle <- readRDS(path_concat(dataDirectory, "lle.RDS"))
kable(lle, row.names = TRUE,
      digits = 3)

## ----tour 2D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center"----
#  p2 <- l_tour(scalingTrainX, color = trainY)
#  l <- l_layer_hull(p2, group = trainY)

## ----2D projection, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=3, fig.height=3, fig.align="center", out.width = "70%"----
include_graphics(path_concat(imageDirectory, "proj2D.PNG"))
proj2D <- readRDS(path_concat(dataDirectory, "proj2D.RDS")) %>% 
  as.matrix()

## ----2D show, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", out.width = "70%"----
#  proj2D <- p2["projection"]

## ----2D projection display, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=3, fig.height=3, fig.align="center", out.width = "70%"----
kable(as.data.frame(proj2D, row.names = colnames(trainX)), 
      digits = 2)

## ----tour 3D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", out.width = "70%"----
#  p3 <- l_tour(scalingTrainX,
#               tour_path = tourr::grand_tour(3),
#               color = trainY,
#               axesLayout = "parallel")
#  proj3D <- p3["projection"]

## ----tour 3D projection, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=3, fig.height=3, fig.align="center", out.width = "70%"----
include_graphics(path_concat(imageDirectory, "proj3D.PNG"))
proj3D <- readRDS(path_concat(dataDirectory, "proj3D.RDS")) %>% 
  as.matrix()

## ----tour 4D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", out.width = "70%"----
#  p4 <- l_tour(scalingTrainX,
#               tour_path = tourr::grand_tour(4),
#               color = trainY,
#               axesLayout = "parallel")
#  proj4D <- p4["projection"]

## ----tour 4D projection, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=3, fig.height=3, fig.align="center", out.width = "70%"----
include_graphics(path_concat(imageDirectory, "proj4D.PNG"))
proj4D <- readRDS(path_concat(dataDirectory, "proj4D.RDS")) %>% 
  as.matrix()

## ----tour 5D, eval = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", out.width = "70%"----
#  p5 <- l_tour(scalingTrainX,
#               tour_path = tourr::grand_tour(5),
#               color = trainY,
#               axesLayout = "parallel")
#  proj5D <- p5["projection"]

## ----tour 5D projection, echo = FALSE, warning=FALSE, message=FALSE, error=FALSE, fig.width=3, fig.height=3, fig.align="center", out.width = "70%"----
include_graphics(path_concat(imageDirectory, "proj5D.PNG"))
proj5D <- readRDS(path_concat(dataDirectory, "proj5D.RDS")) %>% 
  as.matrix()

## ----tour compute, warning=FALSE, message=FALSE, error=FALSE, fig.width=4, fig.height=3, fig.align="center", out.width = "70%"----
tour <- low_dim_knn_pred(list(proj2D, proj3D, 
                              proj4D, proj5D), 
                         fun = function(data, proj) {
                           data %*% as.matrix(proj)
                         },
                         scalingTrainX, trainY,
                         scalingTestX, testY,
                         setNames = FALSE) 
colnames(tour) <- paste0("d = ", 2:5)
kable(tour, row.names = TRUE,
      digits = 3)

## ----bind_states, warning=FALSE, message=FALSE, echo=FALSE, eval = FALSE------
#  callback <- function(W) {
#    # W is the widget path name
#    pred <- knn_pred(scalingTrainX %*% p5["projection"],
#                     trainY,
#                     scalingTestX %*% p5["projection"],
#                     testY,
#                     k = c(5, 10, 20))
#    cat(
#      paste0(
#        "k = 5: accuracy rate ", round(pred[1], 3), "\n",
#        "k = 10: accuracy rate ", round(pred[2], 3), "\n",
#        "k = 20: accuracy rate ", round(pred[3], 3), "\n"
#      )
#    )
#  }
#  l_bind_state(target = l_getPlots(p5),
#               event = "all",
#               callback = callback)
#  # [1] "stateBinding0"

## ----graphical summary, warning=FALSE, message=FALSE, error=FALSE, fig.width=12, fig.height=6, fig.align="center", out.width = "70%"----
rbind(tour, lle, PCA, pChooseD) %>% 
  mutate(k = rep(c(5, 10, 20), 4),
         method = rep(c("tour", "LLE", "PCA", "pChooseD"), each = 3)) %>% 
  pivot_longer(cols = -c(k, method),
               names_to = "Dimensions",
               values_to = "Accuracy") %>%
  mutate(Dimensions = parse_number(Dimensions)) %>%
  ggplot(mapping = aes(x = Dimensions, 
                       y = Accuracy,
                       colour = method)) + 
  geom_path() + 
  facet_wrap(~k) + 
  ggtitle("Facet by the number of neibourhoods")

