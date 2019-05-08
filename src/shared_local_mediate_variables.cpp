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
  
  this->PredictM0 = env["PredictM0"];
  this->PredictM1 = env["PredictM1"];
  this->YModel = env["YModel"];
  
  Rcpp::DataFrame R_y_data = Rcpp::as<Rcpp::DataFrame>(env["y.data"]);
  Rcpp::Language terms = R_y_data.attr("terms");
  Rcpp::CharacterVector term_labels = terms.attr("term.labels");
  
  for(long long int i=0;i<term_labels.size();++i){
    std::string label;
    label = term_labels[i];
    
    if(label == this->treat){
      this->treat_i = i;
    }
    
    if(label == this->mediator){
      this->mediator_i = i;
    }
      
    this->terms.emplace_back(label);
  }
  
  Eigen::Index ncol = R_y_data.cols();
  Eigen::Index nrow = R_y_data.rows();
  
  Eigen::MatrixXd y_data_local(nrow, ncol);
  
  for(long long i=0; i < R_y_data.cols(); ++i ){
    auto col = Rcpp::as<Eigen::VectorXd>(R_y_data[i]);
    y_data_local.col(i) = col;
  }
  
  this->y_data = y_data_local;
  
  this->effects_tmp.reserve(4);
  for(std::size_t i=0;i<4;++i){
    this->effects_tmp.emplace_back(this->n, this->sims);
  }
}

// [[Rcpp::export]]
void test(Rcpp::Environment & env){
  SharedLocalMediateVariables sv;
  sv.initialize_from_environment(env);
  Rcpp::Rcout << sv.y_data.rows() << ',' <<  sv.y_data.cols()<< std::endl;
  Rcpp::Rcout << sv.y_data << std::endl;
}