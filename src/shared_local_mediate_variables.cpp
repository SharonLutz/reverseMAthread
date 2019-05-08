#include "shared_local_mediate_variables.hpp"

// [[Rcpp::plugins(cpp11)]]

SharedLocalMediateVariables::SharedLocalMediateVariables() : 
  PredictM0(0,0),
  PredictM1(0,0),
  YModel(0,0),
  y_data(0,0){}

void SharedLocalMediateVariables::initialize_from_environment(Rcpp::Environment & env) {
  this->tt_switch = {{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}}};
  this->n = env["n"];
  this->sims = env["sims"];
  this->cat_0 = env["cat.0"];
  this->cat_1 = env["cat.1"];
  this->treat = Rcpp::as<std::string>(env["treat"]);
  this->mediator = Rcpp::as<std::string>(env["mediator"]);
  
  Rcpp::NumericMatrix R_PredictM0 = env["PredictM0"];
  Rcpp::NumericMatrix R_PredictM1 = env["PredictM1"];
  Rcpp::NumericMatrix R_YModel = env["YModel"];
  SEXP df_data = env["y.data"];
  Rcpp::DataFrame R_y_data = df_data;
  
  Rcpp::Language terms = R_y_data.attr("terms");
  Rcpp::CharacterVector term_labels = terms.attr("term.labels");
  
  for(long long int i=0;i<term_labels.size();++i){
    std::string label;
    label = term_labels[i];
    this->terms.emplace_back(label);
  }
  
  this->PredictM0 = DoubleMatrix::from_Rmatrix(R_PredictM0);
  this->PredictM1 = DoubleMatrix::from_Rmatrix(R_PredictM1);
  
  this->YModel = DoubleMatrix::from_Rmatrix(R_YModel);
  this->y_data = DoubleMatrix::from_RDataFrame(R_y_data);
  
  this->effects_tmp.reserve(4);
  for(std::size_t i=0;i<4;++i){
    this->effects_tmp.emplace_back(this->n, this->sims);
  }
}