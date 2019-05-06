#include <vector>
#include <array>
#include <cstdint>
#include <sstream>
#include <Rcpp.h>
#include "reverseC_class_defs.hpp"

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

shared_local_mediate_variables shared_local_vars;
void import_vars_from_environment(Rcpp::Environment & env);
// [[Rcpp::export]]
void mediate_helper(Environment &env);

NumericMatrix pred_to_model_mat(DataFrame &pred_mat){
  Rcpp::Language terms = pred_mat.attr("terms");
  CharacterVector term_labels = terms.attr("term.labels");
  GenericVector row_names = pred_mat.attr("row.names");
  CharacterVector names = pred_mat.attr("names");
  
  long long int nrows = row_names.size();
  long long int nterms = term_labels.size();
  
  /*
  for(long long int col_i=0;col_i<nterms;++col_i){
    std::string s = Rcpp::as<std::string>(term_labels[col_i]);
    Rcout << s << std::endl;
  }
  //*/
  
  NumericMatrix model_mat(nrows, nterms+1);
  for(long long int row_i=0;row_i<nrows;++row_i){
    model_mat(row_i,0) = 1.0;
  }
  for(long long int col_i=0;col_i<nterms;++col_i){
    std::string s;
    s = term_labels[col_i];
    NumericVector col = pred_mat[s];
    for(long long int row_i=0;row_i<nrows;++row_i){
      model_mat(row_i,col_i+1) = col[row_i];
    }
  }
  return model_mat;
}

NumericMatrix cpp_matrix_mult(NumericMatrix& m1, NumericMatrix& m2){
  IntegerVector dim1 = m1.attr("dim");
  IntegerVector dim2 = m2.attr("dim");
  unsigned long long cols_1 = dim1[1];
  unsigned long long rows_2 = dim2[0];
  NumericMatrix result(dim1[0], dim2[1]);
  if(cols_1 != rows_2){
    std::stringstream ss;
    ss << "The dimensions of these matrices are not compatible for matrix multiplication: dim1: [";
    ss << dim1[0] << "," << dim1[1] << "]; dim2: [" << dim2[0] << "," << dim2[1] << "]";
    Rf_error(ss.str().c_str());
  }
  for(long long int result_row_i=0;result_row_i<dim1[0];++result_row_i){
    for(long long int result_col_i=0; result_col_i<dim2[1];++result_col_i){
      result(result_row_i, result_col_i) = sum(m1(result_row_i,_) * m2(_,result_col_i));
    }
  }
  
  return result;
}

NumericMatrix cpp_matrix_mult_of_transposed(NumericMatrix& m1, NumericMatrix& m2){
  IntegerVector dim1 = m1.attr("dim");
  IntegerVector dim2 = m2.attr("dim");
  unsigned long long cols_1 = dim1[0];
  unsigned long long rows_2 = dim2[1];
  NumericMatrix result(dim1[1], dim2[0]);
  if(cols_1 != rows_2){
    std::stringstream ss;
    ss << "The dimensions of these matrices are not compatible for matrix multiplication: dim1: [";
    ss << dim1[1] << "," << dim1[0] << "]; dim2: [" << dim2[1] << "," << dim2[0] << "]";
    Rf_error(ss.str().c_str());
  }
  for(long long int result_row_i=0;result_row_i<dim1[1];++result_row_i){
    for(long long int result_col_i=0; result_col_i<dim2[0];++result_col_i){
      result(result_row_i, result_col_i) = sum(m1(_,result_row_i) * m2(result_col_i,_));
    }
  }
  
  return result;
}

void inner_loop(
    std::array<int, 4>& tt,
    std::size_t j, 
    NumericMatrix &Pr1, 
    NumericMatrix &Pr0,
    int cat_t,
    int cat_t_ctrl,
    int cat_c,
    int cat_c_ctrl
    ) {
  Rcpp::Rcout << "inner: " << j << std::endl;
  
  DataFrame pred_data_t = Rcpp::clone(shared_local_vars.y_data);
  
  DataFrame pred_data_c = Rcpp::clone(shared_local_vars.y_data);
  
  pred_data_t[shared_local_vars.treat] = cat_t;
  pred_data_c[shared_local_vars.treat] = cat_c;

  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  NumericVector PredictMt = (shared_local_vars.PredictM1(j,_) * tt[2]) + (shared_local_vars.PredictM0(j,_) * (1 - tt[2]));
  NumericVector PredictMc = (shared_local_vars.PredictM1(j,_) * tt[3]) + (shared_local_vars.PredictM0(j,_) * (1 - tt[3]));
  
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  NumericVector pred_data_t_med_vec = pred_data_t[shared_local_vars.mediator];
  NumericVector pred_data_c_med_vec = pred_data_c[shared_local_vars.mediator];
  
  for(int i=0; i < pred_data_t_med_vec.size(); ++i){
    pred_data_t_med_vec[i] = PredictMt[i];
  }
  for(int i=0; i < pred_data_c_med_vec.size(); ++i){
    pred_data_c_med_vec[i] = PredictMc[i];
  }
  
  NumericMatrix ymat_t = pred_to_model_mat(pred_data_t);
  NumericMatrix ymat_c = pred_to_model_mat(pred_data_c);
  
  NumericVector YModel_j_vec = shared_local_vars.YModel(j,_);

  //NumericMatrix YModel_j_as_matrix(YModel_j_vec.size(),1);
  //YModel_j_as_matrix(_,YModel_j_vec.size()) = YModel_j_vec;
  
  NumericMatrix YModel_j_as_matrix(YModel_j_vec.size(), 1, YModel_j_vec.begin());
  
  Pr1(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_t);
  Pr0(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_c);
  
}

bool compare_vect(const IntegerVector &v1, const IntegerVector &v2){
  if(v1.length() != v2.length()){
    return false;
  }
  for(long long int i=0;i<v1.length();++i){
    if(v1[i] != v2[i]){
      return false;
    }
  }
  return true;
}

NumericMatrix subtract_matrices(const NumericMatrix& mat1, const NumericMatrix &mat2){
  if(!compare_vect(mat1.attr("dim"),mat2.attr("dim"))){
    Rf_error("Matrix Dimensions do not match");
  }
  IntegerVector dimv=mat1.attr("dim");
  NumericMatrix result(dimv[0],dimv[1]);
  long long int total_len = dimv[0]*dimv[1];
  for(long long int i=0;i<total_len;++i){
    result[i] = mat1[i] - mat2[i];
  }
  return result;
}

void outer_loop(std::size_t e) {
  Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = shared_local_vars.tt_switch[e];
  NumericMatrix Pr1(shared_local_vars.n, shared_local_vars.sims);
  NumericMatrix Pr0(shared_local_vars.n, shared_local_vars.sims);
  int cat_t = tt[0] ? shared_local_vars.cat_1 : shared_local_vars.cat_0;
  int cat_t_ctrl = tt[1] ? shared_local_vars.cat_1 : shared_local_vars.cat_0;
  int cat_c = tt[0] ? shared_local_vars.cat_0 : shared_local_vars.cat_1;
  int cat_c_ctrl = tt[1] ? shared_local_vars.cat_0 : shared_local_vars.cat_1;
  
  for(int j=0; j < shared_local_vars.sims; ++j){
    inner_loop(tt, j, Pr1, Pr0, cat_t, cat_t_ctrl, cat_c, cat_c_ctrl);
  }
  
  shared_local_vars.effects_tmp[e] = subtract_matrices(Pr1, Pr0);
  Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
}

void mediate_helper(Environment &env){
  Rcpp::Rcout<< "mediate_helper begin" << std::endl;
  import_vars_from_environment(env);
  for(std::size_t e=0;e<4;++e){
    outer_loop(e);
  }
  env["et1"] = shared_local_vars.effects_tmp[0];
  env["et2"] = shared_local_vars.effects_tmp[1];
  env["et3"] = shared_local_vars.effects_tmp[2];
  env["et4"] = shared_local_vars.effects_tmp[3];
  Rcpp::Rcout<< "mediate_helper end" << std::endl;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


void import_vars_from_environment(Rcpp::Environment &env){
  Rcpp::Rcout<< "import_vars_from_environment begin" << std::endl;
  shared_local_vars.tt_switch = {{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}}};
  shared_local_vars.n = env["n"];
  shared_local_vars.sims = env["sims"];
  shared_local_vars.cat_0 = env["cat.0"];
  shared_local_vars.cat_1 = env["cat.1"];
  shared_local_vars.PredictM0 = Rcpp::as<NumericMatrix>(env["PredictM0"]);
  shared_local_vars.PredictM1 = Rcpp::as<NumericMatrix>(env["PredictM1"]);
  shared_local_vars.YModel = Rcpp::as<NumericMatrix>(env["YModel"]);
  shared_local_vars.y_data = env["y.data"];
  shared_local_vars.treat = Rcpp::as<std::string>(env["treat"]);
  shared_local_vars.mediator = Rcpp::as<std::string>(env["mediator"]);
  for(std::size_t i=0;i<4;++i){
    Rcpp::NumericMatrix new_mat(shared_local_vars.n, shared_local_vars.sims);
    shared_local_vars.effects_tmp.push_back(new_mat);
  }
  Rcpp::Rcout<< "import_vars_from_environment end" << std::endl;
}
