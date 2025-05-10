SharedNeighborDistance = function(Data, k = 5, NThreads = NULL, ComputationInR = FALSE){
  #
  # V = SharedNeighborDistance(Data)
  # V = SharedNeighborDistance(Data, k = 7)
  #
  # DESCRIPTION
  # Calculate the shared neighbor distance for a given dataset.
  #
  # INPUT
  # Data[1:n,1:d]    Numeric matrix with n cases and d dimensions
  # k                Integer defining the number of nearest neighbors
  # ComputationInR    Boolean (Default ComputationInR = FALSE). If FALSE, do
  #                   computation in Rcpp, else in R (very slow).
  #
  # OPTIONAL
  #
  # OUTPUT
  # DistanceMatrix    Numeric matrix with nxn elements
  #
  # DETAILS
  # https://github.com/albert-espin/snn-clustering/blob/master/SNN/snn.py
  #
  # Authors: Quirin Stier, 09/2024
  #

  # Require: parallel, parallelDist, RcppParallel, Rcpp
  x = Sys.time()
  DM         = as.matrix(parallelDist::parDist(Data, method = "euclidean"))
  if(is.null(NThreads)){
    NumThreads = defaultNumThreads() - 1
  }else{
    NumThreads = NThreads
  }
  cl         = makeCluster(NumThreads)
  VtmpRes    = parLapply(cl = cl, X = 1:dim(DM)[1], fun = function(i, localDM, localK){
    order(localDM[i,])[2:(localK+1)]
  }, DM, k)
  NeighborIdx    = do.call(rbind, VtmpRes)
  #NeighborIdx2   = apply(DM, 1, order)[, 1:(k+1)] # Wrong results!
  AdjMat         = matrix(0, nrow = dim(Data)[1], ncol = dim(Data)[1])
  tmpVar         = cbind(rep(1:dim(Data)[1], k), as.vector(NeighborIdx))
  AdjMat[tmpVar] = 1

  if(isTRUE(ComputationInR)){
    DistMat        = matrix(0, nrow = dim(Data)[1], ncol = dim(Data)[1])
    for(i in 2:dim(Data)[1]){
      for(j in 1:(i-1)){
        # Shared neighbor distance
        tmpVar       = 1 - (length(which((AdjMat[j,] + AdjMat[i,]) == 2))/k)
        DistMat[i,j] = tmpVar
        DistMat[j,i] = tmpVar
      }
    }
  }else{
    DistMat = SharedNeighborDistance_Rcpp(AdjMat = AdjMat, N = dim(AdjMat)[1], NumNeighbors = k)
  }
  y = Sys.time()
  Tdiff = y-x
  print(paste0("Computational time: ", round(Tdiff, 2), " ", units(Tdiff)))
  return(DistMat)
}
