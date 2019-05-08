#ifndef STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#define STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP
#include <Rcpp.h>
#include "double_matrix.hpp"
#include <mutex>

// [[Rcpp::plugins(cpp11)]]

struct SharedLocalMediateVariables {
  int n;
  int sims;
  int cat_0;
  int cat_1;
  std::string treat;
  std::string mediator;
  std::vector<std::string> terms;
  std::array< std::array<int, 4>, 4 > tt_switch;
  DoubleMatrix PredictM0;
  DoubleMatrix PredictM1;
  DoubleMatrix YModel;
  DoubleMatrix y_data;
  std::mutex effects_tmp_mutex;
  std::vector<DoubleMatrix>effects_tmp;
  SharedLocalMediateVariables();
  void initialize_from_environment(Rcpp::Environment & env);
};
#endif //STRUCT_SHAREDLOCALMEDIATEVARIABLES_HPP