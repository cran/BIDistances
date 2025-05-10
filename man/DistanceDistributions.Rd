\name{DistanceDistributions}
\alias{DistanceDistributions}
\title{ Distance Distribution}
\description{
Calculates the distribution of the distances between the data points
}
\usage{
DistanceDistributions(Data, DistanceMethods=c('bhjattacharyya', 'bray',
                                                       'canberra', 'chord',
                                                       'divergence', 'euclidean',
                                                       'minkowski', 'geodesic',
                                                       'hellinger', 'kullback',
                                                       'manhattan', 'maximum',
                                                       'soergel', 'wave',
                                                       'whittaker'),
                               CosineNonParallel = TRUE, CorrelationDist = TRUE,
                               Mahalanobis = FALSE, Podani = FALSE,
                               PlotIt = FALSE, PlotSampleSize = 5e3)
}
\arguments{
\item{Data}{[1:n, 1:m] A matrix, containing data as rows.}
\item{DistanceMethods}{Character vector stating all distance methods such as 'euclidean'.}
\item{CosineNonParallel}{Boolean stating if cosine should be computed in parallel.}
\item{CorrelationDist}{Boolean stating if CorrelationDist should be computed.}
\item{Mahalanobis}{Boolean stating if Mahalanobis should be computed.}
\item{Podani}{Boolean stating if Podani should be computed.}
\item{PlotIt}{Boolean: TRUE => create plot. FALSE => no plot.}
\item{PlotSampleSize}{Integer stating the number of samples for plotting.}
}
\value{
List with elements
\item{DistanceMatrix}{[1:n, 1:n] numeric matrix containing the distance matrix}
\item{DistanceChoice}{[1:n, 1:n] numeric matrix containing the distance matrix}
\item{OrderedDistances}{[1:n, 1:n] numeric matrix containing the distance matrix}
\item{ggobject}{ggplot object}
}
\author{Michael Thrun}

\examples{
iris=datasets::iris
if(requireNamespace("DataVisualizations",quietly=TRUE)){
library(DataVisualizations)
DistanceDistributions(as.matrix(iris[,1:4]), c("euclidean"), PlotIt = FALSE)
}
}
