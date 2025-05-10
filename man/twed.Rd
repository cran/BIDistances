\name{twed}
\alias{twed}
\title{twed}
\description{twed}
\usage{
twed(Values1, Values2, Time1, Time2, Nu = 1, Lambda = 1, Degree = 2)
}
\arguments{
\item{Values1}{[1:N1] Numeric vector with values of the first time series.}
\item{Values2}{[1:N1] Numeric vector with values of the second time series.}
\item{Time1}{[1:N1] Numeric vector with time stamps of the first time series.}
\item{Time2}{[1:N1] Numeric vector with time stamps of the second time series.}
\item{Nu}{Optional, Numeric: Elasticity parameter - nu >=0 needed for distance
measure.}
\item{Lambda}{Optional, Numeric: Penalty for deletion operation.}
\item{Degree}{Optional, Integer: Degree of the p norm for local cost.}
}
\value{
List with elements
\item{TWED}{TWED distance between time series Values1 (Time1) and Values2
(Time2)}
\item{DPMatrix}{[1:n, 1:m] Numeric matrix}
}
\author{Quirin Stier}

\examples{
twed(1:10, 1:10, 1:10, 1:10)
}
