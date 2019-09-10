#pragma once

#ifndef STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#define STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include <string>
#include <array>
#include <vector>

struct SharedLocalMediateVariables {
  int n;
  int sims;
  int cat_0;
  int cat_1;
  std::string treat;
  std::string mediator;
  Eigen::Index treat_i;
  Eigen::Index mediator_i;
  std::vector<std::string> terms;
  std::array< std::array<int, 4>, 4 > tt_switch;
  Eigen::MatrixXd PredictM0;
  Eigen::MatrixXd PredictM1;
  Eigen::MatrixXd YModel;
  Eigen::MatrixXd y_data;
  
  Eigen::MatrixXd et1;
  Eigen::MatrixXd et2;
  Eigen::MatrixXd et3;
  Eigen::MatrixXd et4;
  SharedLocalMediateVariables();
  void initialize_from_environment(Rcpp::Environment & env);
  void store_result_diff(Eigen::MatrixXd Pr1, Eigen::MatrixXd Pr0, std::size_t e);
  void export_results(Rcpp::Environment & env);
};
#endif //STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP