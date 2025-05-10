nearest <- function (Data, i, defined) {
  #MT: Korrektur fuer Standardfall
  if(missing(defined)){
    defined = rep(1, ncol(Data))
  }
  distList = Dist2All(Data[i,], Data, defined, knn=2)
  nnind = distList[[2]][2] # nearest neighbor is the point itself, so choose second nearest neighbor
  return(nnind = nnind)
}
