#include "mediate_with_cpp.hpp"
#include "double_matrix.hpp"
#include "shared_local_mediate_variables.hpp"
#include "vector_math_operators.hpp"
#include <thread>
#include <future>
#include <deque>

// [[Rcpp::plugins(cpp11)]]

DoubleMatrix pred_to_model_mat(SharedLocalMediateVariables& sv, DoubleMatrix &pred_mat){
  unsigned long long int nrows = pred_mat.get_num_rows();
  unsigned long long int nterms = sv.terms.size();
  
  DoubleMatrix model_mat(nrows, nterms+1);
  std::vector<std::string> col_names{"(intercept)"};
  
  for(std::size_t i=0;i<nterms;++i){
    col_names.push_back(sv.terms[i]);
  }
  
  model_mat.set_col_names(col_names);
  
  for(unsigned long long int row_i=0;row_i<nrows;++row_i){
    model_mat(row_i,0) = 1.0;
  }
  
  for(unsigned long long int col_i=0;col_i<nterms;++col_i){
    std::string term_label;
    term_label = sv.terms[col_i];
    model_mat.assign_column(pred_mat[term_label], col_i);
  }
  
  return model_mat;
}

bool compare_vect(const Rcpp::IntegerVector &v1, const Rcpp::IntegerVector &v2){
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

//*
void inner_loop(
    SharedLocalMediateVariables& sv,
    std::array<int, 4>& tt,
    std::size_t j, 
    DoubleMatrix &Pr1, 
    DoubleMatrix &Pr0,
    int cat_t,
    int cat_t_ctrl,
    int cat_c,
    int cat_c_ctrl
) {
  // Rcpp::Rcout << "inner: " << j << std::endl;
  
  DoubleMatrix pred_data_t = sv.y_data;
  
  DoubleMatrix pred_data_c = sv.y_data;
  
  // Rcpp::Rcout << "TEST0" << std::endl;
  
  pred_data_t.assign_column(cat_t, sv.treat);
  pred_data_c.assign_column(cat_c, sv.treat);
  
  // Rcpp::Rcout << "TEST1" << std::endl;
  
  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  
  double tt_val_2 = tt[2];
  double tt_val_3 = tt[3];
  
  // Rcpp::Rcout << "sv.PredictM1(j) size: " << sv.PredictM1(j).size() << std::endl;
  // Rcpp::Rcout << "sv.PredictM0(j) size: " << sv.PredictM0(j).size() << std::endl;
  
  std::vector<double> PredictM1_row_j = wrap_DoubleMatrix_column_or_row(sv.PredictM1(j));
  std::vector<double> PredictM0_row_j = wrap_DoubleMatrix_column_or_row(sv.PredictM0(j));
  
  // Rcpp::Rcout << "PredictM1_row_j size: " << PredictM1_row_j.size() << std::endl;
  // Rcpp::Rcout << "PredictM0_row_j size: " << PredictM0_row_j.size() << std::endl;
  
  std::vector<double> PredictMt = (PredictM1_row_j * tt_val_2) + (PredictM0_row_j * (1-tt_val_2));
  std::vector<double> PredictMc = (PredictM1_row_j * tt_val_3) + (PredictM0_row_j * (1-tt_val_3));
  //NumericVector PredictMt = (sv.PredictM1(j,_) * tt[2]) + (sv.PredictM0(j,_) * (1 - tt[2]));
  //NumericVector PredictMc = (sv.PredictM1(j,_) * tt[3]) + (sv.PredictM0(j,_) * (1 - tt[3]));
  
  // Rcpp::Rcout << "PredictMt size: " << PredictMt.size() << std::endl;
  // Rcpp::Rcout << "PredictMc size: " << PredictMc.size() << std::endl;
  
  //Rcpp::Rcout << "TEST2" << std::endl;
  
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  
  pred_data_t.assign_column(PredictMt, sv.mediator);
  pred_data_c.assign_column(PredictMc, sv.mediator);
  
  //Rcpp::Rcout << "TEST3" << std::endl;
  
  DoubleMatrix ymat_t = pred_to_model_mat(sv, pred_data_t);
  DoubleMatrix ymat_c = pred_to_model_mat(sv, pred_data_c);
  
  //Rcpp::Rcout << "TEST4" << std::endl;
  
  std::vector<double> YModel_row_j =  wrap_DoubleMatrix_column_or_row(sv.YModel(j));
  
  DoubleMatrix YModel_j_as_matrix = DoubleMatrix::from_vector(YModel_row_j);
  
  //Rcpp::Rcout << "TEST5" << std::endl;
  
  YModel_j_as_matrix.transpose();
  ymat_t.transpose();
  ymat_c.transpose();
  
  //Rcpp::Rcout << "TEST6" << std::endl;
  
  DoubleMatrix mmult1 = YModel_j_as_matrix * ymat_t;
  DoubleMatrix mmult2 = YModel_j_as_matrix * ymat_t;
  
  //Rcpp::Rcout << "TEST7" << std::endl;
  
  // Rcpp::Rcout << mmult1.get_num_rows() << ',' << mmult1.get_num_cols() << std::endl;
  // Rcpp::Rcout << mmult2.get_num_rows() << ',' << mmult2.get_num_cols() << std::endl;
  
  //Pr1(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_t);
  //Pr0(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_c);
}

void outer_loop(SharedLocalMediateVariables& sv, std::size_t e) {
  // Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = sv.tt_switch[e];
  DoubleMatrix Pr1(sv.n, sv.sims);
  DoubleMatrix Pr0(sv.n, sv.sims);
  int cat_t = tt[0] ? sv.cat_1 : sv.cat_0;
  int cat_t_ctrl = tt[1] ? sv.cat_1 : sv.cat_0;
  int cat_c = tt[0] ? sv.cat_0 : sv.cat_1;
  int cat_c_ctrl = tt[1] ? sv.cat_0 : sv.cat_1;
  
  for(int j=0; j < sv.sims; ++j){
    inner_loop(sv, tt, j, Pr1, Pr0, cat_t, cat_t_ctrl, cat_c, cat_c_ctrl);
  }
  std::lock_guard<std::mutex> lg(sv.effects_tmp_mutex);
  sv.effects_tmp[e] = Pr1 - Pr0;
  // Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
}
//*/

void mediate_helper(Rcpp::Environment &env){
  // Rcpp::Rcout<< "mediate_helper begin" << std::endl;
  //shared_vars.reset(new SharedLocalMediateVariables(env));
  SharedLocalMediateVariables shared_vars;
  shared_vars.initialize_from_environment(env);
  
  //*
  for(std::size_t e=0;e<4;++e){
    outer_loop(shared_vars, e);
  }
  //*/
  //*
  env["et1"] = shared_vars.effects_tmp[0].to_Matrix();
  env["et2"] = shared_vars.effects_tmp[1].to_Matrix();
  env["et3"] = shared_vars.effects_tmp[2].to_Matrix();
  env["et4"] = shared_vars.effects_tmp[3].to_Matrix();
  //*/
  // Rcpp::Rcout<< "mediate_helper end" << std::endl;
  //shared_vars.reset();
}

void threaded_mediate_helper_2_threads(SharedLocalMediateVariables &sv){
  std::vector<std::thread> threads;
  std::vector<std::promise<bool>> promises(4);
  
  auto outer_loop_thread = [&sv, &promises](std::size_t e){
    outer_loop(sv, e);
    promises[e].set_value(true);
  };
  //first 2
  for(std::size_t e=0;e<2;++e){
    threads.emplace_back(outer_loop_thread, e);
  }
  
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  for(std::size_t e=0;e<2;++e){
    auto f = promises[e].get_future();
    f.wait();
  }
  for(std::size_t e=0;e<2;++e){
    threads[e].join();
  }
  
  //second 2
  for(std::size_t e=2;e<4;++e){
    threads.emplace_back(outer_loop_thread, e);
  }
  
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  for(std::size_t e=2;e<4;++e){
    auto f = promises[e].get_future();
    f.wait();
  }
  for(std::size_t e=2;e<4;++e){
    threads[e].join();
  }
  
}

void threaded_mediate_helper_4_threads(SharedLocalMediateVariables &sv){
  std::vector<std::thread> threads;
  std::vector<std::promise<bool>> promises(4);
  
  auto outer_loop_thread = [&sv, &promises](std::size_t e){
    outer_loop(sv, e);
    promises[e].set_value(true);
  };
  
  for(std::size_t e=0;e<4;++e){
    threads.emplace_back(outer_loop_thread, e);
  }
  
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  for(std::size_t e=0;e<4;++e){
    auto f = promises[e].get_future();
    f.wait();
  }
  for(std::size_t e=0;e<4;++e){
    threads[e].join();
  }
}

void outer_loop_with_threads(SharedLocalMediateVariables& sv, std::size_t e, long long int num_threads) {
  std::vector<std::thread> threads;
  std::vector<std::promise<bool>> promises(num_threads);
  std::size_t next_index_value = 0;
  std::mutex remaining_indexes_mutex;
  
  
  
  // Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = sv.tt_switch[e];
  DoubleMatrix Pr1(sv.n, sv.sims);
  DoubleMatrix Pr0(sv.n, sv.sims);
  int cat_t = tt[0] ? sv.cat_1 : sv.cat_0;
  int cat_t_ctrl = tt[1] ? sv.cat_1 : sv.cat_0;
  int cat_c = tt[0] ? sv.cat_0 : sv.cat_1;
  int cat_c_ctrl = tt[1] ? sv.cat_0 : sv.cat_1;
  
  auto inner_loop_thread = [&](std::size_t thread_index){
    std::size_t j;
    while(next_index_value < sv.sims){
      {
        std::lock_guard<std::mutex> lg(remaining_indexes_mutex);
        if(next_index_value >= sv.sims){
          break;
        }
        j=next_index_value++;
      }
      inner_loop(sv, tt, j, Pr1, Pr0, cat_t, cat_t_ctrl, cat_c, cat_c_ctrl);
    }
    promises[thread_index].set_value(true);
  };
  
  for(int ti=0; ti < num_threads; ++ti){
    threads.emplace_back(inner_loop_thread, ti);
  }
  
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  for(std::size_t i=0;i<num_threads;++i){
    auto f = promises[i].get_future();
    f.wait();
  }
  
  for(std::size_t i=0;i<num_threads;++i){
    threads[i].join();
  }
  
  std::lock_guard<std::mutex> lg(sv.effects_tmp_mutex);
  sv.effects_tmp[e] = Pr1 - Pr0;
  // Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
}

void threaded_mediate_helper_N_threads(SharedLocalMediateVariables &sv, long long int num_threads){
  for(std::size_t e=0;e<4;++e){
    outer_loop_with_threads(sv, e, num_threads);
  }
}

void threaded_mediate_helper(Rcpp::Environment &env, long long int num_threads){
  SharedLocalMediateVariables shared_vars;
  shared_vars.initialize_from_environment(env);
  switch(num_threads){
  case 2:
    threaded_mediate_helper_2_threads(shared_vars);
    break;
  case 4:
    threaded_mediate_helper_4_threads(shared_vars);
    break;
  default:
    threaded_mediate_helper_N_threads(shared_vars, num_threads);
    break;
  }
  
  env["et1"] = shared_vars.effects_tmp[0].to_Matrix();
  env["et2"] = shared_vars.effects_tmp[1].to_Matrix();
  env["et3"] = shared_vars.effects_tmp[2].to_Matrix();
  env["et4"] = shared_vars.effects_tmp[3].to_Matrix();
  
  
}

Rcpp::NumericMatrix cpp_mult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2){
  DoubleMatrix dm1 = DoubleMatrix::from_Rmatrix(m1);
  DoubleMatrix dm2 = DoubleMatrix::from_Rmatrix(m2);
  return (dm1 * dm2).to_Matrix();
}

Rcpp::NumericMatrix cpp_tmult(Rcpp::NumericMatrix m1, Rcpp::NumericMatrix m2){
  DoubleMatrix dm1 = DoubleMatrix::from_Rmatrix(m1);
  DoubleMatrix dm2 = DoubleMatrix::from_Rmatrix(m2);
  dm1.transpose();
  dm2.transpose();
  return (dm1 * dm2).to_Matrix();
}
