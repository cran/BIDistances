\name{EuclideanDistances_GPU}
\alias{EuclideanDistances_GPU}
\title{EuclideanDistances_GPU}
\description{EuclideanDistances_GPU}
\usage{
EuclideanDistances_GPU(Data, Weights,OutputType="mat", ctx = NULL)
}
\arguments{
\item{Data}{[1:n,1:d] numeric matrix with data.}

\item{OutputType}{"mat": returns [1:n,1:n] matrix of distances, "dist": returns an object of class "dist"."vec" returns vector of distances}
\item{Weights}{Optional, [1:d] numeric matrix with data.}
\item{ctx}{{Optional,] numeric matrix with data.}
}
}
\value{
Depending on input argument \code{OutputType} either [1:n,1:n] distance matrix or object of class dist or vector of distances
}
\author{Quirin Stier, Michael Thrun}
\examples{
if(isTRUE(require("OpenCL"))){
  if(Sys.info()[1] == "Windows"){
    tryCatch({
      dyn.load("C:/Windows/System32/OpenCL.dll", FALSE, TRUE)
    }, error = function(e) {
      warning("EuclideanDistances_GPU: OpenCL.dll could not be loaded with error:")
      warning(e)
    })
  }
  library(OpenCL)
  AllDevices = unlist(lapply(unlist(OpenCL::oclPlatforms()), OpenCL::oclDevices))
  Device     = AllDevices[length(AllDevices)][[1]]
  ctx        = OpenCL::oclContext(device = Device, precision = "double")
  Data       = as.matrix(iris[,1:4])
  Weights    = rep(1, 4)
  oclDM      = EuclideanDistances_GPU(Data = Data, Weights = Weights, ctx = ctx)
  }
}
\seealso{
\code{\link{as.dist}}
}
\keyword{Distance}
\keyword{Distancematrix}
\keyword{Euclidean}
\keyword{OpenCL}
\keyword{GPU}
\keyword{HPC}
