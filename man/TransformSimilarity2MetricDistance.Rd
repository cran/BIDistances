\name{TransformSimilarity2MetricDistance}
\alias{TransformSimilarity2MetricDistance}
\title{TransformSimilarity2MetricDistance}
\description{TransformSimilarity2MetricDistance}
\usage{
TransformSimilarity2MetricDistance(Similarity)
}
\arguments{
\item{Similarity}{Similarity}
}
\value{
Similarity
}
\author{Michael Thrun}
\examples{
Data_S = fastPdist(as.matrix(iris[,1:4]))
Data_S = Data_S-min(Data_S)
Data_S = Data_S/max(Data_S)
diag(Data_S) = 1
TransformSimilarity2MetricDistance(Data_S)
}
