// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
#include <Rcpp.h>
#include <RcppEigen.h>
#include  "mediation_helper.hpp"

// [[Rcpp::export]]
void mediate_helper(Rcpp::Environment &env){
  // Eigen::setNbThreads(1);
  MediationHelper h(env);
  h();
  // Eigen::setNbThreads(0);
}

// [[Rcpp::export]]
void mediate_helper_variable_exporter(Rcpp::Environment &env){
  // Eigen::setNbThreads(1);
  MediationHelper h(env, true);
  h();
  // Eigen::setNbThreads(0);
}

// [[Rcpp::export]]
void threaded_mediate_helper(Rcpp::Environment &env, long long int num_threads){
  // Eigen::setNbThreads(1);
  if(num_threads > 1){
    Eigen::initParallel();
    MediationHelper h(env, num_threads);
    h();
  } else {
    mediate_helper(env);
  }
  // Eigen::setNbThreads(0);
}

