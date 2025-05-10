## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BIDistances)

## -----------------------------------------------------------------------------
data(Hepta) 
distMatrix = CosinusDistance(Hepta$Data)


## -----------------------------------------------------------------------------
data(Hepta)
V = Dist2All(Hepta$Data[1,],Hepta$Data, method = "euclidean", knn=3)
# Vector of distances from Hepta$Data[1,] to all other rows in Hepta$Data
print(V$distToAll)
# Vector of the indices of the k-nearest-neighbors, according to the euclidean distance
print(V$KNN)

## -----------------------------------------------------------------------------
data(Hepta)
Dmatrix = DistanceMatrix(Hepta$Data, method='euclidean')

## -----------------------------------------------------------------------------
Dmatrix = DistanceMatrix(Hepta$Data, method='minkowski', dim=3)

## -----------------------------------------------------------------------------
data(Hepta)
distMatrix = FractionalDistance(Hepta$Data, p = 1/2)

## -----------------------------------------------------------------------------
data(Hearingloss_N109)
V = Tfidf_dist(Hearingloss_N109$FeatureMatrix_Gene2Term, tf_fun = mean)
# Get distances
dist = V$dist
# Get weights
TfidfWeights = V$TfidfWeights

