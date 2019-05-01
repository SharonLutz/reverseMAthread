#include <vector>
#include <array>

#include <Rcpp.h>

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

// [[Rcpp::export]]
NumericMatrix pred_to_model_mat(DataFrame pred_mat){
  Rcpp::Language terms = pred_mat.attr("terms");
  CharacterVector term_labels = terms.attr("term.labels");
  GenericVector row_names = pred_mat.attr("row.names");
  CharacterVector names = pred_mat.attr("names");
  long long int nrows = row_names.size();
  long long int nterms = term_labels.size();
  
  for(long long int col_i=0;col_i<nterms;++col_i){
    Rcout << term_labels[col_i] << std::endl;
  }
  
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

// [[Rcpp::export]]
NumericMatrix cpp_matrix_transpose(NumericMatrix mat){
  /*
  IntegerVector dim = mat.attr("dim");
  NumericMatrix result(dim[1],dim[0]);
  for(long long int row_i=0;row_i< dim[0];++row_i){
    result(_, row_i) = mat(row_i,_);
  }
  return result;
   */
  return transpose(mat);
}

// [[Rcpp::export]]
NumericMatrix cpp_matrix_mult(NumericMatrix& m1, NumericMatrix& m2){
  IntegerVector dim1 = m1.attr("dim");
  IntegerVector dim2 = m2.attr("dim");
  NumericMatrix result(dim1[0], dim2[1]);
  if(dim1[0] != dim2[1] || dim1[1] != dim2[0]){
    Rf_error("The dimensions of these matrices are not compatible for matrix multiplication");
  }
  for(long long int result_row_i=0;result_row_i<dim1[0];++result_row_i){
    for(long long int result_col_i=0; result_col_i<dim2[1];++result_col_i){
      result(result_row_i, result_col_i) = sum(m1(result_row_i,_) * m2(_,result_col_i));
    }
  }
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix cpp_matrix_mult_of_transposed(NumericMatrix& m1, NumericMatrix& m2){
  IntegerVector dim1 = m1.attr("dim");
  IntegerVector dim2 = m2.attr("dim");
  //tdim1[0] = dim1[1]
  //tdim1[1] = dim1[0]
  //tdim2[0] = dim2[1]
  //tdim2[1] = dim2[0]
  NumericMatrix result(dim1[1], dim2[0]);
  if(dim1[1] != dim2[0] || dim1[0] != dim2[1]){
    Rf_error("The dimensions of these matrices are not compatible for matrix multiplication");
  }
  for(long long int result_row_i=0;result_row_i<dim1[1];++result_row_i){
    for(long long int result_col_i=0; result_col_i<dim2[0];++result_col_i){
      result(result_row_i, result_col_i) = sum(m1(_,result_row_i) * m2(result_col_i,_));
    }
  }
  
  return result;
}

std::vector< std::vector<long long int> > tt_switch{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}};
  
// [[Rcpp::export]]
void inner_loop(Environment env) {
  
  int j = env["j"];
  int cat_t = env["cat.t"];
  int cat_t_ctrl = env["cat.t.ctrl"];
  int cat_c = env["cat.c"];
  int cat_c_ctrl = env["cat.c.ctrl"];
  std::string treat = env["treat"];
  std::string mediator = env["mediator"];
  
  NumericVector tt = env["tt"];
                   
  NumericMatrix Pr1 = env["Pr1"];
  NumericMatrix Pr0 = env["Pr0"];
  
  DataFrame y_data(env["y_data"]);
  
  NumericMatrix PredictM1 = env["PredictM1"];
  NumericMatrix PredictM0 = env["PredictM0"];
  
  NumericMatrix ymat_t = env["ymat.t"];
  NumericMatrix ymat_c = env["ymat.c"];
  
  NumericMatrix ymodel = env["YModel"];
  
  DataFrame pred_data_t = y_data;
  DataFrame pred_data_c = y_data;
  
  pred_data_t[treat] = cat_t;
  pred_data_c[treat] = cat_c;
  
  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  NumericVector PredictMt = (PredictM1(j,_) * tt[2]) + (PredictM0(j,_) * (1 - tt[2]));
  NumericVector PredictMc = (PredictM1(j,_) * tt[3]) + (PredictM0(j,_) * (1 - tt[3]));
  
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  pred_data_t[mediator] = PredictMt;
  pred_data_c[mediator] = PredictMc;
  
  NumericMatrix mmat_t = pred_to_model_mat(pred_data_t);
  NumericMatrix mmat_c = pred_to_model_mat(pred_data_c);
  
  
  NumericVector m1_vec = ymodel(j,_); 
  NumericMatrix m1(m1_vec.size(),1);
  m1(_,m1_vec.size()) = m1_vec;
  
  Pr1(_,j) = cpp_matrix_mult_of_transposed(m1,ymat_t);
  Pr0(_,j) = cpp_matrix_mult_of_transposed(m1,ymat_c);
  
}


// [[Rcpp::export]]
void outer_loop(Environment env) {
  
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


