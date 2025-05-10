// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
struct CalcIVD : public Worker{
  const RVector<double> InputVector;
  RMatrix<double> DistMat;

  CalcIVD(const NumericVector InputVector,
          NumericMatrix DistMat):
    InputVector(InputVector), DistMat(DistMat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = 1; i < end; i++){

      for(std::size_t j = 0; j < i; j++){

        double ValDiff = abs(InputVector[i] - InputVector[j]);

        if(InputVector[i] < InputVector[j]){
          DistMat(i,j) = -ValDiff;
          DistMat(j,i) = ValDiff;
        }else{
          DistMat(i,j) = ValDiff;
          DistMat(j,i) = -ValDiff;
        }

      }
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
NumericMatrix InnerVariableDifferencesRcpp_helper(NumericVector InputVector,
                                                  int N,
                                                  NumericMatrix DistMat){
  CalcIVD distMat(InputVector, DistMat);
  parallelFor(1, N, distMat);
  return DistMat;
}

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
NumericMatrix InnerVariableDifferencesRcpp(NumericVector InputVector, int N){

  NumericMatrix DistMat(N,N);
  DistMat = InnerVariableDifferencesRcpp_helper(InputVector, N, DistMat);

  return(DistMat);
}

