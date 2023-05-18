EuclideanDistances_GPU = function(Data, Weights, OutputType = "mat", ctx = NULL){
  #
  #
  #
  N               = dim(Data)[1]
  DIM             = dim(Data)[2]
  if (missing(Weights)) {
    Weights = rep(1, DIM)
  }

  if(isTRUE(requireNamespace("OpenCL"))) {
    if (Sys.info()[1] == "Windows") {
      tryCatch({
        dyn.load("C:/Windows/System32/OpenCL.dll", FALSE, TRUE)
      }, error = function(e) {
        warning("EuclideanDistances_GPU: OpenCL.dll could not be loaded with error:")
        warning(e)
      })
    }
    if (is.null(ctx)) {
      ctx = OpenCL::oclContext(precision = "best")
    }
    Data = as.matrix(Data)
    if (!is.matrix(Data)) {
      stop("Data must be a matrix.")
    }
    if (!is.vector(Weights)) {
      stop("Weights must be a vector.")
    }

    if (length(Weights) != DIM) {
      stop("Length of vector weights must be equal to the dimension of the data matrix.")
    }

    code <- .ocl_system_file("WeightedEuclideanOCL.cl")
    kernel <- withCallingHandlers({
      OpenCL::oclSimpleKernel(ctx, "weighted_euclidean", code)
    }, warning = function(condition) {
      txt <-
        "OpenCL implementation does not support out-of-order execution"
      if (startsWith(conditionMessage(condition), txt))
        invokeRestart("muffleWarning")
    })

    # Prepare OpenCL kernel
    Size     = N * N
    VecInput = OpenCL::as.clBuffer(as.vector(Data), ctx)
    VecWgts  = OpenCL::as.clBuffer(as.vector(Weights), ctx)
    # Execute OpenCL kernel
    tmpDM    = as.numeric(OpenCL::oclRun(
      kernel = kernel,
      size = Size,
      VecInput,
      VecWgts,
      N,
      DIM,
      dim = c(N, N)
    ))
  }else{
    warning(
      "EuclideanDistances_GPU: OpenCL could not be loaded, falling back to parallel computation on CPU."
    )
    Weights = sqrt(Weights)
    if (!requireNamespace("parallel")) {
      MaxThreads = parallel::detectCores()
      MaxThreads = MaxThreads - 1
      cl = parallel::makeCluster(spec = MaxThreads)
      ScaledData = parallel::parSapply(
        cl = cl,
        X = 1:dim(Data)[2],
        FUN = function(i, Matrix, Weights) {
          Matrix[, i] * Weights[i]
        },
        Data,
        Weights
      )
      tmpDM = parallelDist::parallelDist(x = ScaledData, method = "euclidean")
    } else{
      ScaledData = sapply(
        X = 1:dim(Data)[2],
        FUN = function(i, Matrix, Weights) {
          Matrix[, i] * Weights[i]
        },
        Data,
        Weights
      )
      tmpDM = parallelDist::parallelDist(x = ScaledData, method = "euclidean")
    }
  }
  switch(OutputType,
         "mat" = {
           Distance = matrix(tmpDM, nrow = N, ncol = N)
         },
         "dist" = {
           Distance = as.dist(matrix(tmpDM, nrow = N, ncol = N))
         },
         "vec" = {
           Distance = as.vector(tmpDM)
         },
         {
           Distance = matrix(tmpDM, nrow = N, ncol = N)
           warning("EuclideanDistances_GPU: Incorrect OutputType selected. returning mat.")

         })
  return(Distance)
}
