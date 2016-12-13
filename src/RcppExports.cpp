// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// wedgeRcpp
NumericVector wedgeRcpp(NumericVector a1, NumericVector b1, NumericVector a2, NumericVector b2, int size_vect, double tau, int N, bool lower_tail);
RcppExport SEXP wedge_wedgeRcpp(SEXP a1SEXP, SEXP b1SEXP, SEXP a2SEXP, SEXP b2SEXP, SEXP size_vectSEXP, SEXP tauSEXP, SEXP NSEXP, SEXP lower_tailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< int >::type size_vect(size_vectSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< bool >::type lower_tail(lower_tailSEXP);
    rcpp_result_gen = Rcpp::wrap(wedgeRcpp(a1, b1, a2, b2, size_vect, tau, N, lower_tail));
    return rcpp_result_gen;
END_RCPP
}
