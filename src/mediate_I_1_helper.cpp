#include <Rcpp.h>
#include <vector>
#include <cstdint>
using namespace Rcpp;

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


void coerce_to_factor(IntegerVector &x, CharacterVector levels){
  x.attr("levels") = levels;
  x.attr("class") = "factor";
}


void coerce_to_matrix(NumericVector &x, double rows, double cols) {
  NumericVector d{rows,cols};
  x.attr("dim") = d;
}

std::vector<NumericVector> tt_list{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}};
// [[Rcpp::export]]
void mediation_helper(Environment& env){
  auto n = static_cast<std::size_t>(env["n"]);
  auto sims = static_cast<std::size_t>(env["sims"]);
  DataFrame y_data = env["y.data"];
  DataFrame pred_data_t = env["pred.data.t"];
  DataFrame pred_data_c = env["pred.data.c"];
  SEXP covariates_ref = env["covariates"];
  double cat_0 = env["cat.0"];
  double cat_1 = env["cat.1"];
  bool is_factor_T = env["isFactorT"];
  for(std::size_t e=0; e < 4;++e){
    auto& tt = tt_list[e];
    NumericMatrix Pr1(n, sims);
    NumericMatrix Pr0(n, sims);
    for(std::size_t j=0; j < sims; ++j){
      pred_data_t = pred_data_c = y_data;
      
      if(Rf_isNull(covariates_ref)){
        DataFrame covariates = covariates_ref;
        for(long int p = 0; p < covariates.length(); ++p){
          GenericVector cov_p_as_gv = covariates[p];
          StringVector vl = cov_p_as_gv.names();
          if(Rf_isFactor(pred_data_t[_,vl])){
            StringVector levels = as<IntegerVector>(y_data[_,vl]).attr("levels");
            IntegerVector cov_p_as_iv=as<IntegerVector>(cov_p_as_gv);
            coerce_to_factor(cov_p_as_iv, levels);
          }else{
            pred_data_t[_,vl] = pred_data_c[_,vl] = cov_p_as_gv;
          }
        }
      }
      
      double cat_t;
      double cat_t_ctrl;
      double cat_c;
      double cat_c_ctrl;
      if(tt[0]){
        cat_t = cat_1;
        cat_t_ctrl = cat_0;
      } else {
        cat_t = cat_0;
        cat_t_ctrl = cat_1;
      }
      if(tt[1]){
        cat_c = cat_1;
        cat_c_ctrl = cat_0;
      } else {
        cat_c = cat_0;
        cat_c_ctrl = cat_1;
      }
      
      if(is_factor_T){
        StringVector t_levels = env["t.levels"];
        if(Rf_isNull(env["control"])){
          
        }
      } else {
        
      }
      
    }
  }
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
