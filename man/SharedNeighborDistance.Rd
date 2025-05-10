\name{SharedNeighborDistance}
\alias{SharedNeighborDistance}
\title{Shared Neighbor Distance}
\description{Calculates the Shared Neighbor Distance}
\usage{SharedNeighborDistance(Data, k = 5, NThreads = NULL, ComputationInR = FALSE)}
\arguments{
\item{Data}{[1:n,1:d] matrix with n cases, d variables}
\item{k}{Integer defining the number of nearest neighbors}
\item{NThreads}{Number of threads in parallel computation.}
\item{ComputationInR}{Boolean (Default ComputationInR = FALSE). If FALSE, do
computation in Rcpp, else in R (very slow).}
}
\value{
\item{Distance}{[1:n,1:n] symmetric matrix, containing the distanes of the cases (rows) for the given data}
}
\references{
\url{https://github.com/albert-espin/snn-clustering/blob/master/SNN/snn.py}
}
\author{Quirin Stier}
\examples{
data(Hepta)
distMatrix = SharedNeighborDistance(Hepta$Data, NThreads = 1, ComputationInR=TRUE)
}

\keyword{Distance}
\concept{Shared Neighbor Distance}
\concept{Bioinformatic Distance}
