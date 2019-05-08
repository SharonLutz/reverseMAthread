#ifndef MEDIATE_WITH_CPP_HPP
#define MEDIATE_WITH_CPP_HPP
#include "double_matrix.hpp"
#include "shared_local_mediate_variables.hpp"
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
void mediate_helper(Rcpp::Environment &env);

void threaded_mediate_helper(Rcpp::Environment &env);

DoubleMatrix pred_to_model_mat(SharedLocalMediateVariables& sv, DoubleMatrix &pred_mat);
Rcpp::NumericMatrix cpp_mult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);
Rcpp::NumericMatrix cpp_tmult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);

#endif //MEDIATE_WITH_CPP_HPP