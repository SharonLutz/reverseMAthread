#ifndef MEDIATE_WITH_CPP_HPP
#define MEDIATE_WITH_CPP_HPP
#include "shared_local_mediate_variables.hpp"
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]

void mediate_helper(Rcpp::Environment &env);

void threaded_mediate_helper(Rcpp::Environment &env, long long int num_threads);

Eigen::MatrixXd pred_to_model_mat(SharedLocalMediateVariables& sv, Eigen::MatrixXd &pred_mat);


#endif //MEDIATE_WITH_CPP_HPP