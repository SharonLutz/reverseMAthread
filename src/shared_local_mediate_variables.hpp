#ifndef STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#define STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#include <Rcpp.h>
#include <RcppEigen.h>
#include <mutex>

// [[Rcpp::plugins(cpp11)]]

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
  std::mutex effects_tmp_mutex;
  std::vector<Eigen::MatrixXd> effects_tmp;
  SharedLocalMediateVariables();
  void initialize_from_environment(Rcpp::Environment & env);
};
#endif //STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP