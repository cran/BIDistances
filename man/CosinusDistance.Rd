\name{CosinusDistance}
\alias{CosinusDistance}
\title{Cosine Distance}
\description{Calculates the cosine distance}
\usage{CosinusDistance(Data)}
\arguments{
\item{Data}{
[1:n,1:d] matrix with n cases, d variables
}
}
\details{
\url{https://en.wikipedia.org/wiki/Cosine_similarity}
}
\value{
\item{Distance}{[1:n,1:n] symmetric matrix, containing the distanes of the cases (rows) for the given data}
}
\author{
Michael Thrun
}
\note{
The cosine distance is calculated by calculating the cosine similarity \eqn{d(i,j)=\max{s}-s(i,j)}, where \eqn{s} is the cosine similarity and the \eqn{d} the cosine distance.
}

\examples{
data(Hepta)
distMatrix = CosinusDistance(Hepta$Data)
}

\keyword{Distance}
\keyword{Cosine}
\concept{Cosine Distance}
\concept{Bioinformatic Distance}
