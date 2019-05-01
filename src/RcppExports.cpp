// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mediation_helper
void mediation_helper(Environment& env);
RcppExport SEXP _reverseC_mediation_helper(SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment& >::type env(envSEXP);
    mediation_helper(env);
    return R_NilValue;
END_RCPP
}
// cpp_matrix_transpose
NumericMatrix cpp_matrix_transpose(NumericMatrix& mat);
RcppExport SEXP _reverseC_cpp_matrix_transpose(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_matrix_transpose(mat));
    return rcpp_result_gen;
END_RCPP
}
// cpp_matrix_mult
NumericMatrix cpp_matrix_mult(NumericMatrix& m1, NumericMatrix& m2);
RcppExport SEXP _reverseC_cpp_matrix_mult(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_matrix_mult(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// cpp_matrix_mult_of_transposed
NumericMatrix cpp_matrix_mult_of_transposed(NumericMatrix& m1, NumericMatrix& m2);
RcppExport SEXP _reverseC_cpp_matrix_mult_of_transposed(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_matrix_mult_of_transposed(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// inner_loop_before_model_matrix
void inner_loop_before_model_matrix(Environment env);
RcppExport SEXP _reverseC_inner_loop_before_model_matrix(SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    inner_loop_before_model_matrix(env);
    return R_NilValue;
END_RCPP
}
// inner_loop_after_model_matrix
void inner_loop_after_model_matrix(Environment env);
RcppExport SEXP _reverseC_inner_loop_after_model_matrix(SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    inner_loop_after_model_matrix(env);
    return R_NilValue;
END_RCPP
}
// outer_loop_before_inner_loop
void outer_loop_before_inner_loop(Environment env);
RcppExport SEXP _reverseC_outer_loop_before_inner_loop(SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    outer_loop_before_inner_loop(env);
    return R_NilValue;
END_RCPP
}
// outer_loop_after_inner_loop
void outer_loop_after_inner_loop(Environment env);
RcppExport SEXP _reverseC_outer_loop_after_inner_loop(SEXP envSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type env(envSEXP);
    outer_loop_after_inner_loop(env);
    return R_NilValue;
END_RCPP
}
// derp
NumericVector derp(NumericVector x);
RcppExport SEXP _reverseC_derp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(derp(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_reverseC_mediation_helper", (DL_FUNC) &_reverseC_mediation_helper, 1},
    {"_reverseC_cpp_matrix_transpose", (DL_FUNC) &_reverseC_cpp_matrix_transpose, 1},
    {"_reverseC_cpp_matrix_mult", (DL_FUNC) &_reverseC_cpp_matrix_mult, 2},
    {"_reverseC_cpp_matrix_mult_of_transposed", (DL_FUNC) &_reverseC_cpp_matrix_mult_of_transposed, 2},
    {"_reverseC_inner_loop_before_model_matrix", (DL_FUNC) &_reverseC_inner_loop_before_model_matrix, 1},
    {"_reverseC_inner_loop_after_model_matrix", (DL_FUNC) &_reverseC_inner_loop_after_model_matrix, 1},
    {"_reverseC_outer_loop_before_inner_loop", (DL_FUNC) &_reverseC_outer_loop_before_inner_loop, 1},
    {"_reverseC_outer_loop_after_inner_loop", (DL_FUNC) &_reverseC_outer_loop_after_inner_loop, 1},
    {"_reverseC_derp", (DL_FUNC) &_reverseC_derp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_reverseC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
