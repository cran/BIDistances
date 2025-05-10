\name{DistanceMatrix}
\alias{DistanceMatrix}
\alias{pdist}
\title{Pairwise distance between pairs of objects}
\description{
computes the distance between objects in the data matrix, X, using the method specified by method
}
\usage{
DistanceMatrix(X,method='euclidean',dim=2,outputisvector=FALSE)
}
\arguments{
  \item{X}{data matrix [1:n,1:d], n cases d variables}
  \item{method}{Optional, method specified by distance string:
                'binary','canberra','cityblock','euclidean,
                'sqEuclidean',
								'maximum','cosine','chebychev','jaccard,'kendallM','kendallD'
								'mahalanobis','minkowski','manhattan','braycur','cosine','wasserstein','pearsonD','spearmanD','pearsonM','spearmanM'
	}
  \item{dim}{Optional: if method="minkowski", or wasserstein, choose scalar. For minkowski the ppth root of the sum of the ppth powers of the differences of the components. For wasserstein the order, default should be then 1}
  \item{outputisvector}{Optional: should the output be converted to a vector}
}
\details{
If possible uses implementation parallelized by the parallelDist package. Otherwise R implementations besides Euclidean for which a GPU implementation is provided.

'binary'
(aka asymmetric binary): The vectors are regarded as binary bits,
 so non-zero elements are 'on' and zero elements are 'off'.
 The distance is the proportion of bits in which only one is
 on amongst those in which at least one is on.


'cityblock'==manhattan

'maximum': Maximum distance between two components of x and y (supremum norm)

'cosine'
calculates a similarity matrix sim between all column vectors of a matrix x.
This matrix might be a document-term matrix, so columns would be expected
 to be documents and rows to be terms.
the distances is than defined with D=max(sim)-sim

'jaccard'
Jaccard index is computed as 2B/(1+B), where B is Bray-Curtis dissimilarity:
the number of items which occur in both elements divided by the total
number of items in the elements (Sneath, 1957). This measure is often
also called: binary, asymmetric binary, etc.

'mahalanobis'
 the squared generalized Mahalanobis distance between all pairs of
 rows in a data frame with respect to a covariance matrix.
 The element of the i-th row and j-th column of the distance matrix is defined as
\eqn{D_{ij}^2 = (\bold{x}_i - \bold{x}_j)' \bold{S}^{-1} (\bold{x}_i - \bold{x}_j)}

'minkowski':The p norm, the pth root of the sum of the pth powers
 of the differences of the components.

'manhattan': Absolute distance between the two vectors (1 norm aka L_1).

'chebychev'=max(abs(x-y)),

'canberra'=sum abs(x-y)/sum(abs(x)-abs(y)),
Terms with zero numerator and denominator are omitted from the sum
and treated as if the values were missing.
This is intended for non-negative values (e.g., counts): taking
the absolute value of the denominator is a 1998 R modification to
avoid negative distances.

'braycur'=sum abs(x -y)/abs(x+y)

'pearsonM' metric, see [Legendre, 1986] or [Bock,1974, pp.77-79]
sqrt((1 - r)+1)/2) with r beeing the Pearson's correlation coefficient.

'spearmanM' metric, see [Legendre, 1986] or [Bock,1974, pp.77-79]
sqrt((1 - r)+1)/2)  with r beeing Spearman's correlation coefficient.

'kendallM' metric, see [Legendre, 1986] or [Bock,1974, pp.77-79]
sqrt((1 - r)+1)/2) with tau beeing Kendalls's correlation coefficient.

'pearsonD' dissimilarity
1 - r with r beeing the Pearson's correlation coefficient.

'spearmanD' dissimilarity
1 - r with r beeing Spearman's correlation coefficient.

'kendallD' dissimilarity
1 - r with tau beeing Kendalls's correlation coefficient.


'cosine'
s. wiki for similarity conversion: max(S)-S(i,j)
}
\value{
	\item{Dmatrix}{[1:n,1:n] Distance Marix: Pairwise distance
	between pairs of objects}
}
\references{
Sneath, P. H. A. (1957) Some thoughts on bacterial classification. Journal of General Microbiology 17, pages 184-200.

Leydesdorff, L. (2005) Similarity Measures, Author Cocitation Analysis,and Information Theory. In: JASIST 56(7), pp.769-772.

Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.

Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979) Multivariate Analysis. Academic Press.

Borg, I. and Groenen, P. (1997) Modern Multidimensional Scaling. Theory and Applications. Springer.

Mahalanobis, P. C. (1936) On the generalized distance in statistics. Proceedings of The National Institute of Sciences of India, 12:49-55.
}
\author{Michael Thrun}
\examples{
  data(Hepta)
	Dmatrix = DistanceMatrix(Hepta$Data,method='euclidean')
}
