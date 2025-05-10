\name{WassersteinDist}
\alias{WassersteinDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wasserstein Distance
}
\description{
Computes the Wasserstein distance for a data matrix
}
\usage{
WassersteinDist(Data, p = 1, InverseWeighting = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
data matrix of n cases and d feautures
}
  \item{p}{
scalar higher than one, the power to which the Euclidean distance between points is taken in order to compute transportation costs.
}
  \item{InverseWeighting}{
weighting per row can be either 1 (FALSE) or 1/n (TRUE)
}
}
\details{
Wasserstein distance, also known as Earth Mover’s Distance (EMD) is the distance between two probability distributions over a region D. The Wasserstein distance of order p is defined as the p-th root of the total cost incurred when transporting measure a to measure b in an optimal way, where the cost of transporting a unit of mass from
x to  y is given as the p-th power of the Euclidean distance.

It is claimed to be useful for distributions that do not align well with traditional measures like Euclidean distance.
}
\value{
matrix of distances, symmetric
}
\references{
...
}
\author{
Michae Thrun
}



\seealso{
\code{\link[transport]{wasserstein1d}}
}
\examples{
data(Hepta)
distMat=WassersteinDist(Hepta$Data)
}
\keyword{wasserstein}

\concept{ wasserstein distance }

