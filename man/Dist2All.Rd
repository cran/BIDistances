\name{Dist2All}
\alias{Dist2All}
\title{ Distances to all data points}
\description{
Calculates all distances from a given vector to the rows of a matrix.
}
\usage{
Dist2All(X, Data, SelectFeatures, method = "euclidean",p=2,knn=1)
}
\arguments{
\item{X}{ A vector containing the data point to be compared to data.}
\item{Data}{ A matrix containing the data points to be compared with x.}
\item{SelectFeatures}{ A vector of the same length as x and the rows of data, containing TRUE for all columns of the data to be compared and any other value for columns to be discarded.}
\item{method}{
(Optional) String marking, which distance measure is to be used. Euclidean by default.
}
\item{p}{
(Optional) Scalar, The pp-th root of the sum of the pp-th powers of the differences of the components. Default is 2
}
\item{knn}{
(Optional) Scalar, gives the number of the indices of the k nearest neighbors returned. Default is 1
}
}
\value{
List with
\item{distToAll}{A vector containing the distances from x to all rows of data.}
\item{KNN}{Numeric vector, containing the indices of the k nearest neighbors (rows) to the given points}
}
\note{This function is very inefficient for large Data.}
\author{Michael Thrun}

\examples{
data(Hepta)
Dist2All(Hepta$Data[1,],Hepta$Data)
}

\keyword{Distance}
\concept{Bioinformatic Distance}
\concept{Nearest Neighbors}

