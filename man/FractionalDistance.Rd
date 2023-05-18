\name{FractionalDistance}
\alias{FractionalDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates fractional distances
}
\description{
Calculates distance matrix, through \eqn{\left( \sum_{i=1}^{n} |x_i - y_i|^p \right)^{1/p}}
}
\usage{
FractionalDistance(Data, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
[1:n,1:d] Matrix, with n cases, d variables
}
  \item{p}{
Scalar, value for p
}
}
\details{
Values of p < 1 can be used, which can be useful for high-dimensional data, see references.
}
\value{
\item{DistanceMatrix}{[1:n,1:n] symmetric Matrix, containing the distances between the cases (rows) of the input matrix}
}
\references{
Aggrawal, C. C., Hinneburg, A., Keim, D. (2001), On the Suprising Behavior of Distance Metrics in High Dimensional Space.
}
\author{
Michael Thrun
}

\examples{
data(Hepta)
distMatrix = FractionalDistance(Hepta$Data, p = 1/2)
}
\keyword{Distance}
\keyword{Metric}
\concept{Fractional Distance}
\concept{Bioinformatic Distance}
