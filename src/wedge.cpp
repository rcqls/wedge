#include <Rcpp.h>
#include <math.h>
#include "../inst/include/wedgeC.h"

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector wedgeRcpp(NumericVector a1, NumericVector b1, NumericVector a2, NumericVector b2, int size_vect, double tau, int N, bool lower_tail) {
  NumericVector wp(size_vect);
  for(int i=0;i<size_vect;i++) wp[i]=wedgeC(a1[i],b1[i],a2[i],b2[i],tau,N,lower_tail);
  return wp;
}
