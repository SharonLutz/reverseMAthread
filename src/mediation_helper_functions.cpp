#include <vector>
#include <array>
#include <map>
#include <cstdint>
#include <sstream>
#include <memory>
#include <Rcpp.h>
#include <RcppEigen.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void mediate_helper(Rcpp::Environment &env);

// [[Rcpp::export]]
void threaded_mediate_helper(Rcpp::Environment &env, long long int num_threads);

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_mult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_tmult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2);

// [[Rcpp::export]]
Eigen::VectorXd getEigenValues(Eigen::Map<Eigen::MatrixXd> M) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  return es.eigenvalues();
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_mat_mult(Eigen::Map<Eigen::MatrixXd> M1, Eigen::Map<Eigen::MatrixXd> M2){
  return M1 * M2;
}
