#pragma once

#ifndef REVERSEC_MEDIATION_HELPER_HPP
#define REVERSEC_MEDIATION_HELPER_HPP
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include "shared_local_mediate_variables.hpp"
#include <thread>
#include <future>
#include <deque>
#include <array>
#include <cstdint>
#include <atomic>
#include <mutex>

class OuterLoopVars {
public:
  std::size_t e;
  std::array<int, 4>& tt;
  Eigen::MatrixXd Pr1;
  Eigen::MatrixXd Pr0;
  int cat_t;
  int cat_t_ctrl;
  int cat_c;
  int cat_c_ctrl;
  
  OuterLoopVars(SharedLocalMediateVariables &sv, std::size_t e);
};

class MediationHelper {
protected:
  
  std::size_t num_threads;
  
  Rcpp::Environment &env;
  
  bool export_loop_vars;
  
  bool vars_exported;
  
  SharedLocalMediateVariables sv;
  
  std::vector<std::thread> threads;
  std::vector<std::promise<void>> promises;
  
  void pred_to_model_mat(Eigen::MatrixXd &pred_mat, Eigen::MatrixXd &model_mat);
  
  void inner_loop(OuterLoopVars &olv, std::size_t j );
  void outer_loop(std::size_t e);
  void outer_loop_with_threaded_inner_loop(std::size_t e);
  
  void launcher_n_1();
  void launcher_n_2();
  void launcher_n_4();
  void launcher_n_other();
  
public:
  explicit MediationHelper(Rcpp::Environment &env, bool export_loop_vars=false);
  explicit MediationHelper(Rcpp::Environment &env, long long int num_threads);
  void operator()();
};
#endif //REVERSEC_MEDIATION_HELPER_HPP