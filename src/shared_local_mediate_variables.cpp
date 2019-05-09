// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
#include "shared_local_mediate_variables.hpp"

SharedLocalMediateVariables::SharedLocalMediateVariables() : 
  PredictM0(0,0),
  PredictM1(0,0),
  YModel(0,0),
  y_data(0,0),
  et1(0,0),
  et2(0,0),
  et3(0,0),
  et4(0,0){}

void SharedLocalMediateVariables::initialize_from_environment(Rcpp::Environment & env) {
  this->tt_switch = {{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}}};
  this->n = env["n"];
  this->sims = env["sims"];
  this->cat_0 = env["cat.0"];
  this->cat_1 = env["cat.1"];
  this->treat = Rcpp::as<std::string>(env["treat"]);
  this->mediator = Rcpp::as<std::string>(env["mediator"]);
  
  this->PredictM0 = env["PredictM0"];
  this->PredictM1 = env["PredictM1"];
  this->YModel = env["YModel"];
  
  Rcpp::DataFrame R_y_data = Rcpp::as<Rcpp::DataFrame>(env["y.data"]);
  // Rcpp::Language terms = R_y_data.attr("terms");
  // Rcpp::CharacterVector variable_labels = terms.attr("term.labels");
  
  Rcpp::CharacterVector variable_labels = R_y_data.names();
  
  Eigen::Index ncol = R_y_data.cols();
  Eigen::Index nrow = R_y_data.rows();
  
  Eigen::MatrixXd y_data_local(nrow, ncol);
  
  for(long long int i=0;i<variable_labels.size();++i){
    std::string label;
    label = variable_labels[i];
    
    if(label == this->treat){
      this->treat_i = i;
    }
    if(label == this->mediator){
      this->mediator_i = i;
    }
    long long int j = R_y_data.findName(label);
    
    auto col = Rcpp::as<Eigen::VectorXd>(R_y_data[j]);
    y_data_local.col(i) = col;
    this->terms.emplace_back(label);
  }
  
  this->y_data = y_data_local;
}

// [[Rcpp::export]]
void test(Rcpp::Environment & env){
  SharedLocalMediateVariables sv;
  sv.initialize_from_environment(env);
  Rcpp::Rcout << sv.y_data.rows() << ',' <<  sv.y_data.cols()<< std::endl;
  Rcpp::Rcout << sv.y_data << std::endl;
}

void SharedLocalMediateVariables::store_result_diff(Eigen::MatrixXd Pr1, Eigen::MatrixXd Pr0, std::size_t e){
  // Rcpp::Rcout << (Pr1(0,0) - Pr0(0,0)) << std::endl;
  // Rcpp::Rcout << (Pr1 - Pr0) << std::endl;
  switch(e){
  default:
    Rf_error("invalid value of e");
  case 0:
    et1 = Pr1 - Pr0;
    return;
  case 1:
    et2 = Pr1 - Pr0;
    return;
  case 2:
    et3 = Pr1 - Pr0;
    return;
  case 3:
    et4 = Pr1 - Pr0;
    return;
  }
}

void SharedLocalMediateVariables::export_results(Rcpp::Environment & env){
  env["et1"] = this->et1;
  env["et2"] = this->et2;
  env["et3"] = this->et3;
  env["et4"] = this->et4;
}