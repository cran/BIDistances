\name{ToroidDist2All}
\alias{ToroidDist2All}
\title{Calculate toroid Euclidean Distances }
\arguments{
\item{positionxy}{One datapoint}
\item{AllPositions(1:AnzData:2)}{All Other dataPoints}
\item{Lines,Columns}{Size of planar grid}
}
\value{Dist2All(1:AnzData,1:AnzData);      distance(s) between XY and AllPositions}
\description{
 Calculate toroid Euclidean Distances }
\author{
MT}
\examples{
positionxy = c(1,1)
AllPositions = rbind(c(2,3), c(5,2))
Lines = 40
Columns = 80
ToroidDist2All(positionxy, AllPositions, Lines, Columns)
}
