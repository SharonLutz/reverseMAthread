
#ifndef REVERSEC_CLASS_DEF_HEADER_H
#define REVERSEC_CLASS_DEF_HEADER_H

#include <vector>
#include <array>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

struct shared_local_mediate_variables {
  int n;
  int sims;
  int cat_0;
  int cat_1;
  std::string treat;
  std::string mediator;
  Rcpp::NumericMatrix PredictM0;
  Rcpp::NumericMatrix PredictM1;
  Rcpp::NumericMatrix YModel;
  Rcpp::DataFrame y_data;
  std::array< std::array<int, 4>, 4 > tt_switch;
  Rcpp::NumericMatrix et1;
  Rcpp::NumericMatrix et2;
  Rcpp::NumericMatrix et3;
  Rcpp::NumericMatrix et4;
    
  std::vector<Rcpp::NumericMatrix *>effects_tmp;
};

#endif //REVERSEC_CLASS_DEF_HEADER_H