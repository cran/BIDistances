// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace RcppParallel;
using namespace Rcpp;

using namespace std;

// [[Rcpp::depends(RcppParallel)]]
struct CalcSND : public Worker{
  const RMatrix<double> AdjMat;
  int N;
  double NumNeighbors;
  RMatrix<double> DistMat;

  CalcSND(const NumericMatrix AdjMat,
          int N, double NumNeighbors,
          NumericMatrix DistMat):
    AdjMat(AdjMat), N(N), NumNeighbors(NumNeighbors), DistMat(DistMat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){

      for(std::size_t j = 0; j < i; j++){

        int Counter = 0;

        for(int k = 0; k < N; k++){
          if((AdjMat(i,k) == 1) && (AdjMat(j,k) == 1)){
            Counter = Counter + 1;
          }
        }

        //double tmpVar = Counter;///NumNeighbors;
        DistMat(i,j) = 1-(Counter/NumNeighbors);
        DistMat(j,i) = 1-(Counter/NumNeighbors);
      }
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
NumericMatrix SharedNeighborDistance_helper(NumericMatrix AdjMat, NumericMatrix DistMat, int N, double NumNeighbors){
  CalcSND distMat(AdjMat, N, NumNeighbors, DistMat);
  parallelFor(1, N, distMat);
  return DistMat;
}

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
NumericMatrix SharedNeighborDistance_Rcpp(NumericMatrix AdjMat, int N, double NumNeighbors){

  NumericMatrix DistMat(N,N);
  DistMat = SharedNeighborDistance_helper(AdjMat, DistMat, N, NumNeighbors);

  return(DistMat);
}
