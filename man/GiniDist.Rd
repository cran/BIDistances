\name{GiniDist}
\alias{GiniDist}
\title{GiniDist}
\description{
Calculates pairwise gini distances
}
\usage{
GiniDist(Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Data}{[1:n,1:d] data to calculate distances to
}
}
\value{
dist[1:n,1:n]   distances
}
\author{
Michael Thrun
}
\examples{
GiniDist(as.matrix(iris[,1:4]))
}
