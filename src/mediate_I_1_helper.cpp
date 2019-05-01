#include <vector>
#include <cstdint>
#include <memory>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//[[Rcpp::export]]
void mediation_helper(Environment &env);

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

struct shared_reference_vars{
  std::size_t n;
  std::size_t sims;
  int cat_0;
  int cat_1;
  DataFrame y_data;
  std::vector<NumericMatrix> effects_tmp;
  CharacterVector treat;
  CharacterVector mediator;
  long long int treat_i;
  long long int mediator_i;
  
  shared_reference_vars(Environment &env): 
    n(static_cast<std::size_t>(env["n"])),
    sims(static_cast<std::size_t>(env["sims"])),
    cat_0(env["cat.0"]),
    cat_1(env["cat.1"]),
    y_data(env["y.data"]),
    effects_tmp(),
    treat(env["treat"]),
    mediator(env["mediator"]){
    for(std::size_t i=0;i<4;++i){
      effects_tmp.emplace_back(n,sims);
    }
    CharacterVector ydat_names = y_data.attr("names");
    bool found_treat=false;
    bool found_med=false;
    for(long long int i=0;i<ydat_names.length();++i){
      if(ydat_names[i] == treat[0]){
        found_treat=true;
        treat_i = i;
      }
      if(ydat_names[i] == mediator[0]){
        found_med=true;
        mediator_i = i;
      }
    }
    if(!found_med){
      Rf_error("Did not find treatment index");
    }
    if(!found_treat){
      Rf_error("Did not find mediator index");
    }
  }
};

void assign_to_all(NumericVector& obj, double val){
  for(long long int i=0;i<obj.length();++i){
    obj[i]=val;
  }
}

std::unique_ptr<shared_reference_vars> edv_ptr;

void inner_loop(
    std::size_t j, 
    NumericMatrix &Pr1, 
    NumericMatrix &Pr0,
    int cat_t,
    int cat_t_ctrl,
    int cat_c,
    int cat_c_ctrl){
  DataFrame pred_data_t = edv_ptr->y_data;
  DataFrame pred_data_c = edv_ptr->y_data;
  //pred.data.t[,treat] <- cat.t
  //pred.data.c[,treat] <- cat.c
  //assign_to_all(pred_data_t(_,edv_ptr->treat_i), cat_t);
  //assign_to_all(pred_data_c(_,edv_ptr->treat_i), cat_t);
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

std::vector<std::vector<int>> tt_list{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}};

void outer_loop(std::size_t e){
  std::vector<int>& tt = tt_list[e];
  NumericMatrix Pr1(edv_ptr->n, edv_ptr->sims);
  NumericMatrix Pr0(edv_ptr->n, edv_ptr->sims);
  int cat_t = tt[0] ? edv_ptr->cat_1 : edv_ptr->cat_0;
  int cat_t_ctrl = tt[1] ? edv_ptr->cat_1 : edv_ptr->cat_0;
  int cat_c = tt[0] ? edv_ptr->cat_0 : edv_ptr->cat_1;
  int cat_c_ctrl = tt[1] ? edv_ptr->cat_0 : edv_ptr->cat_1;
  
  for(std::size_t j=0; j < edv_ptr->sims; ++j){
    inner_loop(j, Pr1, Pr0, cat_t, cat_t_ctrl, cat_c, cat_c_ctrl);
  }
  
  edv_ptr->effects_tmp[e] = subtract_matrices(Pr1, Pr0);
}

void check_value(const std::vector<std::size_t>& dims, std::size_t index_value, std::size_t index_pos){
  if(index_pos > dims.size()){
    Rf_error("index position has no equivalent");
  }
  if(index_value >= dims[index_pos]){
    Rf_error("index value too large for dimension");
  }
}

std::size_t get_offset_of(const std::vector<std::size_t>& dims, const std::vector<std::size_t>& indices){
  if(indices.size() != dims.size()){
    Rf_error("wrong number of indices given for this offset");
  }
  check_value(dims, indices[0], 0);
  std::size_t result = indices[0];
  std::size_t factor = 1;
  for(std::size_t i=1;i<dims.size();++i){
    check_value(dims, indices[i], i);
    std::size_t index_value = indices[i];
    factor *= dims[i-1];
    result += index_value * factor;
  }
  return result;
}

NumericMatrix extract_matrix(NumericVector &arr, long long int last_index){
  if(Rf_isNull(arr.attr("dim"))){
    Rf_error("vector has no dim attribute");
  }
  IntegerVector dim_v = arr.attr("dim");
  if(last_index >= dim_v[2]){
    Rf_error("index out of bounds");
  }
  NumericMatrix result(dim_v[0], dim_v[1]);
  std::vector<long long int> indices_of_requested_matrix{0,0,last_index};
  long long int count_elements_in_matrix = dim_v[0] * dim_v[1];
  long long int offset = count_elements_in_matrix * last_index;
  for(long long int i=offset,j=0; i < offset + count_elements_in_matrix ; ++i,++j){
    result[j]=arr[i];
  }
   
  return result;
  
}



void mediation_helper(Environment &env){
  //variables that get used in the loop that are previously defined
  //y.data
  //variables defined in the loop that remain in existence after the loop
  edv_ptr.reset(new shared_reference_vars(env));
  
  for(std::size_t e=0; e < 4;++e){
    outer_loop(e);
  }
  //*/
  env["et1"] = edv_ptr->effects_tmp[0];
  env["et2"] = edv_ptr->effects_tmp[1];
  env["et3"] = edv_ptr->effects_tmp[2];
  env["et4"] = edv_ptr->effects_tmp[3];
  edv_ptr.reset();
}



