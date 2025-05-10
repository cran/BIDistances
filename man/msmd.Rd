\name{msmd}
\alias{msmd}
\title{msmd}
\description{msmd}
\usage{
msmd(Values1, Values2, ParameterC)
}
\arguments{
\item{Values1}{[1:N1] Numeric vector with values of the first time series.}
\item{Values2}{[1:N1] Numeric vector with values of the second time series.}
\item{ParameterC}{Numeric vector with time stamps of the first time series.}
}
\value{
List with elements
\item{Value}{Distance measure}
}
\author{Quirin Stier}
\examples{
msmd(1:10, 1:10)
}

