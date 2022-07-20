// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// MatrixToEntropy
NumericVector MatrixToEntropy(IntegerMatrix x);
RcppExport SEXP _scCDC_MatrixToEntropy(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(MatrixToEntropy(x));
    return rcpp_result_gen;
END_RCPP
}
// VectorToEntropy
double VectorToEntropy(IntegerVector x);
RcppExport SEXP _scCDC_VectorToEntropy(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(VectorToEntropy(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scCDC_MatrixToEntropy", (DL_FUNC) &_scCDC_MatrixToEntropy, 1},
    {"_scCDC_VectorToEntropy", (DL_FUNC) &_scCDC_VectorToEntropy, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_scCDC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}