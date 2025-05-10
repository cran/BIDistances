\name{VariablePrecision}
\alias{VariablePrecision}
\title{VariablePrecision}
\description{Computes the variable precision}
\usage{
VariablePrecision(Variable)
}
\arguments{
\item{Variable}{Numeric vector [1:n] or matrix [1:n, 1:d]}
}
\value{
MinAbsDiff, MinAbsNZDiff, MinExpo
}
\author{
Michael Thrun
}
\examples{
data(Hepta)
distMat = VariablePrecision(as.matrix(iris[, 1]))

distMat = VariablePrecision(as.matrix(iris[, 1:4]))
}
