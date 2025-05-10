\name{Tfidf_dist}
\alias{Tfidf_dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Term frequency-inverse document frequency distance}
\description{
Computes the term frequency inverse document frequency (tfidf) distance for a FeatureMatrix_Gene2GoTerm. In case of genes with annotated GOterms from gene ontology genes can be interpreted as documents and GOterms as terms.
}
\usage{
Tfidf_dist(FeatureMatrix_Gene2GoTerm, tf_fun = mean)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{FeatureMatrix_Gene2GoTerm}{
[1:n,1:d] Matrix, with n genes and d GO-Terms.
}
  \item{tf_fun}{
Function, defining the numerator value in the normalized Term-frequency. The default is the mean of the not 0 values.
}
}
\details{
For the FeatureMatrix_Gene2GoTerm it is: \cr
FeatureMatrix_Gene2GoTerm[i,j] > 0 iff GOterm j is relevant for gene i. The FeatureMatrix_Gene2GoTerm[i,j] > 1 if the specific gene is annotated by in a specific GO-Term with more than one evidence code 	FeatureMatrix_Gene2GoTerm[i,j] is the frequency of term js occurance in document i. \cr
}
\references{
Stier, Q. and Thrun, M., C.: Deriving homogeneous subsets from gene sets by exploiting the Gene Ontology, Informatica, in review, 2023
}

\value{
List with
\item{dist}{Numeric vector containing the tdfidf distances between the documents = absolute difference of TfidfWeights}
\item{TfidfWeights}{[1:n] Numeric vector containing the term frequence inverse document frequency weights used for the distance, given as the Term frequency*Inverse document frequency}
}

\author{Michael Thrun}
\examples{
data(Hearingloss_N109)
V = Tfidf_dist(Hearingloss_N109$FeatureMatrix_Gene2Term)
dist = V$dist
TfidfWeights = V$TfidfWeights
}
\keyword{Distance}
\keyword{TermFrequence}
\concept{Bioinformatic Distance}
\concept{Term frequence}
\concept{Inverse frequency}

