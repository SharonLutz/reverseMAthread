#include "mediate_with_cpp.hpp"
#include "shared_local_mediate_variables.hpp"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <thread>
#include <future>
#include <deque>


// [[Rcpp::plugins(cpp11)]]

void pred_to_model_mat(SharedLocalMediateVariables& sv, Eigen::MatrixXd &pred_mat, Eigen::MatrixXd &model_mat){
  model_mat.col(0).setOnes();
  for(std::size_t col_i=0;col_i<sv.terms.size();++col_i){
    model_mat.col(col_i+1) = pred_mat.col(col_i);
  }
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
    Eigen::MatrixXd &Pr1, 
    Eigen::MatrixXd &Pr0,
    int cat_t,
    int cat_t_ctrl,
    int cat_c,
    int cat_c_ctrl
) {
  // Rcpp::Rcout << "inner: " << j << std::endl;
  
  Eigen::MatrixXd pred_data_t = sv.y_data;
  Eigen::MatrixXd pred_data_c = sv.y_data;
  
  pred_data_t.col(sv.treat_i).setConstant(cat_t);
  pred_data_c.col(sv.treat_i).setConstant(cat_c);
  
  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  
  pred_data_t.col(sv.mediator_i) = (sv.PredictM1.row(j) * tt[2]) + (sv.PredictM0.row(j) * (1-tt[2]));
  pred_data_c.col(sv.mediator_i) = (sv.PredictM1.row(j) * tt[3]) + (sv.PredictM0.row(j) * (1-tt[3]));
  
  Eigen::MatrixXd ymat_t(pred_data_t.rows(), sv.terms.size() + 1);
  Eigen::MatrixXd ymat_c(pred_data_c.rows(), sv.terms.size() + 1);
  
  pred_to_model_mat(sv, pred_data_t, ymat_t);
  pred_to_model_mat(sv, pred_data_c, ymat_c);
  
  //Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
  //Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
  Pr1.col(j) = (sv.YModel.row(j) * ymat_t.transpose()).row(0);
  Pr0.col(j) = (sv.YModel.row(j) * ymat_c.transpose()).row(0);
  
  // Rcpp::Rcout << "TEST3" << std::endl;
}

void outer_loop(SharedLocalMediateVariables& sv, std::size_t e) {
  // Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = sv.tt_switch[e];
  Eigen::MatrixXd Pr1(sv.n, sv.sims);
  Eigen::MatrixXd Pr0(sv.n, sv.sims);
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
  env["et1"] = shared_vars.effects_tmp[0];
  env["et2"] = shared_vars.effects_tmp[1];
  env["et3"] = shared_vars.effects_tmp[2];
  env["et4"] = shared_vars.effects_tmp[3];
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
  std::size_t num_sims = sv.sims;
  
  
  // Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = sv.tt_switch[e];
  Eigen::MatrixXd Pr1(sv.n, sv.sims);
  Eigen::MatrixXd Pr0(sv.n, sv.sims);
  int cat_t = tt[0] ? sv.cat_1 : sv.cat_0;
  int cat_t_ctrl = tt[1] ? sv.cat_1 : sv.cat_0;
  int cat_c = tt[0] ? sv.cat_0 : sv.cat_1;
  int cat_c_ctrl = tt[1] ? sv.cat_0 : sv.cat_1;
  
  auto inner_loop_thread = [&](std::size_t thread_index){
    std::size_t j;
    while(next_index_value < num_sims){
      {
        std::lock_guard<std::mutex> lg(remaining_indexes_mutex);
        if(next_index_value >= num_sims){
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
  
  for(long long int i=0;i<num_threads;++i){
    auto f = promises[i].get_future();
    f.wait();
  }
  
  for(long long int i=0;i<num_threads;++i){
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
  
  env["et1"] = shared_vars.effects_tmp[0];
  env["et2"] = shared_vars.effects_tmp[1];
  env["et3"] = shared_vars.effects_tmp[2];
  env["et4"] = shared_vars.effects_tmp[3];
  
  
}
