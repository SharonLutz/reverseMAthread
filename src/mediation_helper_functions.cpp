#include <vector>
#include <array>
#include <map>
#include <cstdint>
#include <sstream>
#include <memory>
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
void mediate_helper(Environment &env);

class DoubleMatrix{
protected:
  bool col_names_assigned;
  bool row_names_assigned;
  bool transposed;
  std::size_t num_rows;
  std::size_t num_cols;
  std::vector<double> data;
  std::vector<std::vector<double *>> access_via_rows;
  std::vector<std::vector<double *>> access_via_cols;
  std::vector<std::string> row_names;
  std::vector<std::string> col_names;
  std::map<std::string, std::size_t> row_name_to_index;
  std::map<std::string, std::size_t> col_name_to_index;
  
  void initialize_access();
  
public:
  DoubleMatrix(std::size_t num_rows, std::size_t num_cols);
  DoubleMatrix(Rcpp::IntegerVector dims);
  DoubleMatrix(const DoubleMatrix &other);
  DoubleMatrix &operator=(const DoubleMatrix &other);
  
  std::size_t get_size();
  void transpose();
  std::size_t get_num_rows();
  std::size_t get_num_cols();
  
  void set_col_names(Rcpp::CharacterVector col_names);
  void set_row_names(Rcpp::CharacterVector row_names);
  
  void set_col_names(std::vector<std::string> &col_names);
  void set_row_names(std::vector<std::string> &row_names);
  
  std::string get_row_name(std::size_t row_i);
  std::string get_col_name(std::size_t col_i);
  
  double& operator()(std::size_t row_i, std::size_t col_i);
  std::vector<double*> &operator[](std::size_t col_i);
  std::vector<double*> &operator()(std::size_t row_i);
  std::vector<double*> &operator[](std::string &col_name);
  std::vector<double*> &operator()(std::string &row_name);
  std::vector<double*> &operator[](Rcpp::CharacterVector col_name);
  std::vector<double*> &operator()(Rcpp::CharacterVector row_name);
  
  void assign_column(std::vector<double*> col, std::size_t col_i);
  void assign_row(std::vector<double*> row, std::size_t row_i);
  
  void assign_column(std::vector<double*> col, std::string& col_name);
  void assign_row(std::vector<double*> row, std::string& row_name);
  
  DoubleMatrix operator*(DoubleMatrix &other);
  DoubleMatrix operator+(DoubleMatrix &other);
  DoubleMatrix operator-(DoubleMatrix &other);
  DoubleMatrix operator*(double scalar);
  DoubleMatrix operator+(double scalar);
  DoubleMatrix operator-(double scalar);
  DoubleMatrix operator/(double scalar);
  
  static DoubleMatrix from_Rmatrix(Rcpp::NumericMatrix mat);
  static DoubleMatrix from_Rvector(Rcpp::NumericVector vec);
  static DoubleMatrix from_vector(std::vector<double> &svec);
  static DoubleMatrix from_Rvector(Rcpp::NumericVector vec, std::size_t num_rows, std::size_t num_cols);
  static DoubleMatrix from_vector(std::vector<double> &svec, std::size_t num_rows, std::size_t num_cols);
  static DoubleMatrix from_RDataFrame(Rcpp::DataFrame df);
  
  Rcpp::DataFrame to_DataFrame();
  Rcpp::NumericMatrix to_Matrix();
};

DoubleMatrix::DoubleMatrix(std::size_t num_rows, std::size_t num_cols) : 
  col_names_assigned(false),
  row_names_assigned(false),
  transposed(false),
  num_rows(num_rows),
  num_cols(num_cols),
  data(num_rows*num_cols,0){
  initialize_access();
}
DoubleMatrix::DoubleMatrix(Rcpp::IntegerVector dims) : 
  col_names_assigned(false),
  row_names_assigned(false),
  transposed(false),
  num_rows(dims[0]),
  num_cols(dims[1]),
  data(dims[0]*dims[1],0){
  initialize_access();
  }
DoubleMatrix::DoubleMatrix(const DoubleMatrix &other) : 
  col_names_assigned(other.col_names_assigned),
  row_names_assigned(other.row_names_assigned),
  transposed(other.transposed),
  num_rows(other.num_rows),
  num_cols(other.num_cols),
  data(other.data),
  row_names(other.row_names),
  col_names(other.col_names),
  row_name_to_index(other.row_name_to_index),
  col_name_to_index(other.col_name_to_index){
  initialize_access();
  }
DoubleMatrix& DoubleMatrix::operator=(const DoubleMatrix &other){
  this->col_names_assigned = other.col_names_assigned;
  this->row_names_assigned = other.row_names_assigned;
  this->transposed = other.transposed;
  this->num_rows = other.num_rows;
  this->num_cols = other.num_cols;
  this->data = other.data;
  this->row_names = other.row_names;
  this->col_names = other.col_names;
  this->row_name_to_index = other.row_name_to_index;
  this->col_name_to_index = other.col_name_to_index;
  this->initialize_access();
  return (*this);
}

void DoubleMatrix::initialize_access(){
  access_via_rows.resize(this->num_rows);
  access_via_cols.resize(this->num_cols);
  for(std::size_t row_i=0; row_i< this->num_rows; ++row_i){
    auto &row = access_via_rows[row_i];
    row.resize(this->num_cols, nullptr);
    for(std::size_t col_i=0; col_i< this->num_cols; ++col_i){
      row[col_i] = &data[row_i*(this->num_rows) + col_i];
    }
  }
  for(std::size_t col_i=0; col_i < this->num_cols; ++col_i){
    auto &col = access_via_cols[col_i];
    col.resize(this->num_rows, nullptr);
    for(std::size_t row_i=0; row_i < this->num_rows; ++row_i){
      col[row_i] = &data[row_i*this->num_rows + col_i];
    }
  }
}


std::size_t DoubleMatrix::get_size(){
  return this->num_cols * this->num_rows;
}

std::size_t DoubleMatrix::get_num_rows(){
  return this->transposed ? this->num_cols : this->num_rows;
}

std::size_t DoubleMatrix::get_num_cols(){
  return this->transposed ? this->num_rows : this->num_cols;
}

void DoubleMatrix::transpose(){
  this->transposed = !(this->transposed);
}

double& DoubleMatrix::operator()(std::size_t row_i, std::size_t col_i){
  std::size_t row_count = this->get_num_rows();
  std::size_t col_count = this->get_num_cols();
  bool row_oob = row_i >= row_count;
  bool col_oob = col_i >= col_count;
  if(row_oob || col_oob){
    std::stringstream ss;
    ss << "index out of bounds";
    if(row_oob){
      ss << "; rows: " << row_i << " >= " << row_count;
    }
    if(col_oob){
      ss << "; cols: " << col_i << " > = " << col_count;
    }
    Rf_error(ss.str().c_str());
  }
  return this->transposed ? data[col_i * col_count + row_i] : data[row_i * row_count + col_i];
}


std::vector<double*>& DoubleMatrix::operator[](std::size_t col_i){
  return this->transposed ? access_via_rows[col_i] : access_via_cols[col_i];
}

std::vector<double*>& DoubleMatrix::operator()(std::size_t row_i){
  return this->transposed ? access_via_cols[row_i] : access_via_rows[row_i];
}

void DoubleMatrix::set_col_names(Rcpp::CharacterVector col_names){
  std::vector<std::string> names;
  names.reserve(col_names.size());
  for(long long i=0; i < col_names.size(); ++i){
    names.push_back(Rcpp::as<std::string>(col_names[i]));
  }
  this->set_col_names(names);
}

void DoubleMatrix::set_row_names(Rcpp::CharacterVector row_names){
  std::vector<std::string> names;
  names.reserve(row_names.size());
  for(long long i=0; i < row_names.size(); ++i){
    names.push_back(Rcpp::as<std::string>(row_names[i]));
  }
  this->set_row_names(names);
}

void DoubleMatrix::set_col_names(std::vector<std::string> &new_col_names){
  std::size_t col_count = this->get_num_cols();
  if(col_count != new_col_names.size()){
    std::stringstream ss;
    ss << "names should be the same length as the column dimension: " << new_col_names.size() << " != " << col_count;
    Rf_error(ss.str().c_str());
  }
  if(this->transposed){
    this->row_names.clear();
    this->row_names.reserve(col_count);
    this->row_name_to_index.clear();
    for(std::size_t i=0; i < new_col_names.size(); ++i){
      this->row_names.push_back(new_col_names[i]);
    }
  } else {
    this->col_names.clear();
    this->col_names.reserve(col_count);
    this->col_name_to_index.clear();
    for(std::size_t i=0; i < new_col_names.size(); ++i){
      this->col_names.push_back(new_col_names[i]);
    }
  }
}

void DoubleMatrix::set_row_names(std::vector<std::string> &new_row_names){
  std::size_t row_count = this->get_num_rows();
  if(row_count != new_row_names.size()){
    std::stringstream ss;
    ss << "names should be the same length as the row dimension: " << new_row_names.size() << " != " << row_count;
    Rf_error(ss.str().c_str());
  }
  if(this->transposed){
    this->col_names.clear();
    this->col_names.reserve(row_count);
    this->col_name_to_index.clear();
    for(std::size_t i=0; i < new_row_names.size(); ++i){
      this->col_names.push_back(new_row_names[i]);
    }
  } else {
    this->row_names.clear();
    this->row_names.reserve(row_count);
    this->row_name_to_index.clear();
    for(std::size_t i=0; i < new_row_names.size(); ++i){
      this->col_names.push_back(new_row_names[i]);
    }
  }
}

std::string DoubleMatrix::get_row_name(std::size_t row_i){
  if(this->transposed){
    if(this->col_names_assigned){
      return this->col_names[row_i];
    } else {
      std::stringstream ss;
      ss << row_i;
      return ss.str();
    }
  } else {
    if(this->row_names_assigned){
      return this->row_names[row_i];
    } else {
      std::stringstream ss;
      ss << row_i;
      return ss.str();
    }
  }
}

std::string DoubleMatrix::get_col_name(std::size_t col_i){
  if(this->transposed){
    if(this->row_names_assigned){
      return this->row_names[col_i];
    } else {
      std::stringstream ss;
      ss << "V" << col_i;
      return ss.str();
    }
  } else {
    if(this->col_names_assigned){
      return this->col_names[col_i];
    } else {
      std::stringstream ss;
      ss << "V" << col_i;
      return ss.str();
    }
  }
}

std::vector<double*>& DoubleMatrix::operator[](std::string &col_name){
  std::size_t index;
  if(this->transposed){
    if(this->row_name_to_index.count(col_name) == 0){
      
    }
    index = this->row_name_to_index[col_name];
  } else {
    if(this->col_name_to_index.count(col_name) == 0){
      
    }
    index = this->col_name_to_index[col_name];
  }
  return this->operator[](index);
}

std::vector<double*>& DoubleMatrix::operator()(std::string &row_name){
  std::size_t index;
  if(this->transposed){
    if(this->col_name_to_index.count(row_name) == 0){
      
    }
    index = this->col_name_to_index[row_name];
  } else {
    if(this->row_name_to_index.count(row_name) == 0){
      
    }
    index = this->row_name_to_index[row_name];
  }
  return this->operator()(index);
}

std::vector<double*>& DoubleMatrix::operator[](Rcpp::CharacterVector col_name){
  if(col_name.size() != 1){
    std::stringstream ss;
    ss << "too many strings given, expected 1, got: "<< col_name.size();
    Rf_error(ss.str().c_str());
  }
  std::string name = Rcpp::as<std::string>(col_name[0]);
  return this->operator[](name);
}

std::vector<double*>& DoubleMatrix::operator()(Rcpp::CharacterVector row_name){
  if(row_name.size() != 1){
    std::stringstream ss;
    ss << "too many strings given, expected 1, got: "<< row_name.size();
    Rf_error(ss.str().c_str());
  }
  std::string name= Rcpp::as<std::string>(row_name[0]);
  return this->operator()(name);
}

DoubleMatrix DoubleMatrix::operator*(DoubleMatrix &other){
  unsigned long long cols_1 = this->get_num_cols();
  unsigned long long rows_1 = this->get_num_rows();
  unsigned long long cols_2 = other.get_num_cols();
  unsigned long long rows_2 = other.get_num_rows();
  
  if(cols_1 != rows_2){
    std::stringstream ss;
    ss << "The dimensions of these matrices are not compatible for matrix multiplication: dim1: [";
    ss << rows_1 << "," << cols_1 << "]; dim2: [" << rows_2 << "," << cols_2 << "]";
    Rf_error(ss.str().c_str());
  }
  
  DoubleMatrix result(rows_1, cols_2);
  
  for(unsigned long long int result_row_i=0;result_row_i<rows_1;++result_row_i){
    auto &row_1_i = this->operator()(result_row_i);
    for(unsigned long long int result_col_i=0; result_col_i<cols_2;++result_col_i){
      auto &col_2_i = other.operator()(result_col_i);
      double sum_of_vector_mult=0;
      for(std::size_t i=0; i<row_1_i.size(); ++i){
        sum_of_vector_mult += (*row_1_i[i]) * (*col_2_i[i]);
      }
      result(result_row_i, result_col_i) = sum_of_vector_mult;
    }
  }
  
  return result;
}

DoubleMatrix DoubleMatrix::operator+(DoubleMatrix &other){
  bool rows_mismatch = this->get_num_rows() == other.get_num_rows();
  bool cols_mismatch = this->get_num_cols() == other.get_num_cols();
  if(rows_mismatch || cols_mismatch){
    std::stringstream ss;
    ss << "index mismatch between matrices: ";
    ss << '[' << this->get_num_rows() << ',' << this->get_num_cols() <<']';
    ss << " vs ";
    ss << '[' << other.get_num_rows() << ',' << other.get_num_cols() <<']';
    Rf_error(ss.str().c_str());
  }
  DoubleMatrix result(*this);
  for(std::size_t row_i=0; row_i< this->get_num_rows(); ++row_i){
    for(std::size_t col_i=0; col_i < this->get_num_cols(); ++col_i){
      result(row_i, col_i) += other(row_i,col_i);
    }
  }
  return result;
}

DoubleMatrix DoubleMatrix::operator-(DoubleMatrix &other){
  bool rows_mismatch = this->get_num_rows() == other.get_num_rows();
  bool cols_mismatch = this->get_num_cols() == other.get_num_cols();
  if(rows_mismatch || cols_mismatch){
    std::stringstream ss;
    ss << "index mismatch between matrices: ";
    ss << '[' << this->get_num_rows() << ',' << this->get_num_cols() <<']';
    ss << " vs ";
    ss << '[' << other.get_num_rows() << ',' << other.get_num_cols() <<']';
    Rf_error(ss.str().c_str());
  }
  DoubleMatrix result(*this);
  for(std::size_t row_i=0; row_i< this->get_num_rows(); ++row_i){
    for(std::size_t col_i=0; col_i < this->get_num_cols(); ++col_i){
      result(row_i, col_i) -= other(row_i,col_i);
    }
  }
  return result;
}

DoubleMatrix DoubleMatrix::operator*(double scalar){
  DoubleMatrix result(*this);
  for(auto & item : result.data){
    item *= scalar;
  }
  return result;
}

DoubleMatrix DoubleMatrix::operator+(double scalar){
  DoubleMatrix result(*this);
  for(auto & item : result.data){
    item += scalar;
  }
  return result;
}

DoubleMatrix DoubleMatrix::operator-(double scalar){
  DoubleMatrix result(*this);
  for(auto & item : result.data){
    item -= scalar;
  }
  return result;
}

DoubleMatrix DoubleMatrix::operator/(double scalar){
  DoubleMatrix result(*this);
  for(auto & item : result.data){
    item /= scalar;
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_Rmatrix(Rcpp::NumericMatrix mat){
  IntegerVector dim = mat.attr("dim");
  DoubleMatrix result(dim[0], dim[1]);
  for(long long row_i=0; row_i < dim[0];++row_i){
    for(long long col_i=0; col_i < dim[1];++col_i){
      result(row_i, col_i) = mat(row_i, col_i);
    }
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_Rvector(Rcpp::NumericVector vec){
  std::size_t num_rows = vec.length();
  DoubleMatrix result(num_rows, 1);
  for(unsigned long long row_i=0; row_i < num_rows;++row_i){
    result(row_i, 1) = vec[row_i];
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_vector(std::vector<double> &svec){
  std::size_t num_rows = svec.size();
  DoubleMatrix result(num_rows, 1);
  for(unsigned long long row_i=0; row_i < num_rows;++row_i){
    result(row_i, 1) = svec[row_i];
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_Rvector(Rcpp::NumericVector vec, std::size_t num_rows, std::size_t num_cols){
  std::size_t veclen = vec.length();
  if(num_rows * num_cols != veclen){
    std::stringstream ss;
    Rf_error(ss.str().c_str());
  }
  DoubleMatrix result(num_rows, num_cols);
  for(std::size_t row_i=0; row_i < num_rows;++row_i){
    for(std::size_t col_i=0; col_i < num_cols;++col_i){
      result(row_i,col_i) = vec[row_i*num_rows + col_i];
    }
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_vector(std::vector<double> &svec, std::size_t num_rows, std::size_t num_cols){
  if(num_rows * num_cols != svec.size()){
    std::stringstream ss;
    Rf_error(ss.str().c_str());
  }
  DoubleMatrix result(num_rows, num_cols);
  for(std::size_t row_i=0; row_i < num_rows;++row_i){
    for(std::size_t col_i=0; col_i < num_cols;++col_i){
      result(row_i,col_i) = svec[row_i*num_rows + col_i];
    }
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_RDataFrame(Rcpp::DataFrame df){
  IntegerVector dim = df.attr("dim");
  DoubleMatrix result(dim[0], dim[1]);
  CharacterVector col_names = df.attr("names");
  CharacterVector row_names = df.attr("row.names");
  for(long long col_i=0; col_i < dim[1];++col_i){
    Rcpp::NumericVector row_vec = df[col_i];
    for(long long row_i=0; row_i < dim[0];++row_i){
      result(row_i,col_i) = row_vec[row_i];
    }
  }
  result.set_col_names(col_names);
  result.set_row_names(row_names);
  return result;
}

Rcpp::DataFrame DoubleMatrix::to_DataFrame(){
  Rcpp::DataFrame result;
  std::size_t n_cols = this->get_num_cols();
  std::size_t n_rows = this->get_num_rows();
  CharacterVector col_names(n_cols);
  for(std::size_t col_i =0; col_i < n_cols; ++col_i){
    col_names[col_i] = this->get_col_name(col_i);
  }
  for(std::size_t col_i =0; col_i < n_cols; ++col_i){
    Rcpp::NumericVector col_vec(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      col_vec[row_i] = this->operator()(row_i, col_i);
    }
    result[Rcpp::as<CharacterVector>(col_names[col_i])]=col_vec;
  }
  if(this->row_names_assigned){
    CharacterVector row_names(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      row_names[row_i] = this->get_row_name(row_i);
    }
    result.attr("row.names") = row_names;
  } else {
    IntegerVector row_names(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      row_names[row_i]=row_i;
    }
    result.attr("row.names") = row_names;
  }
  return result;
}

Rcpp::NumericMatrix DoubleMatrix::to_Matrix(){
  Rcpp::NumericMatrix result(this->get_num_rows(), this->get_num_cols());
  for(std::size_t row_i=0;row_i < this->get_num_rows(); ++row_i){
    for(std::size_t col_i=0;col_i < this->get_num_cols(); ++col_i){
      result(row_i,col_i) = this->operator()(row_i, col_i);
    }
  }
  return result;
}

struct shared_local_mediate_variables {
  int n;
  int sims;
  int cat_0;
  int cat_1;
  std::string treat;
  std::string mediator;
  DoubleMatrix PredictM0;
  DoubleMatrix PredictM1;
  DoubleMatrix YModel;
  DoubleMatrix y_data;
  std::array< std::array<int, 4>, 4 > tt_switch;
  
  std::vector<DoubleMatrix>effects_tmp;
  shared_local_mediate_variables(Rcpp::Environment & env);
};

std::unique_ptr<shared_local_mediate_variables> shared_vars;

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
/*
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
  
  DataFrame pred_data_t = Rcpp::clone(shared_vars->y_data);
  
  DataFrame pred_data_c = Rcpp::clone(shared_vars->y_data);
  
  pred_data_t[shared_vars->treat] = cat_t;
  pred_data_c[shared_vars->treat] = cat_c;

  //PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
  //PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
  NumericVector PredictMt = (shared_vars->PredictM1(j,_) * tt[2]) + (shared_vars->PredictM0(j,_) * (1 - tt[2]));
  NumericVector PredictMc = (shared_vars->PredictM1(j,_) * tt[3]) + (shared_vars->PredictM0(j,_) * (1 - tt[3]));
  
  //pred.data.t[,mediator] <- PredictMt
  //pred.data.c[,mediator] <- PredictMc
  NumericVector pred_data_t_med_vec = pred_data_t[shared_vars->mediator];
  NumericVector pred_data_c_med_vec = pred_data_c[shared_vars->mediator];
  
  for(int i=0; i < pred_data_t_med_vec.size(); ++i){
    pred_data_t_med_vec[i] = PredictMt[i];
  }
  for(int i=0; i < pred_data_c_med_vec.size(); ++i){
    pred_data_c_med_vec[i] = PredictMc[i];
  }
  
  NumericMatrix ymat_t = pred_to_model_mat(pred_data_t);
  NumericMatrix ymat_c = pred_to_model_mat(pred_data_c);
  
  NumericVector YModel_j_vec = shared_vars->YModel(j,_);

  //NumericMatrix YModel_j_as_matrix(YModel_j_vec.size(),1);
  //YModel_j_as_matrix(_,YModel_j_vec.size()) = YModel_j_vec;
  
  NumericMatrix YModel_j_as_matrix(YModel_j_vec.size(), 1, YModel_j_vec.begin());
  
  Pr1(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_t);
  Pr0(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_c);
  
}

void outer_loop(std::size_t e) {
  Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
  std::array<int, 4>& tt = shared_vars->tt_switch[e];
  NumericMatrix Pr1(shared_vars->n, shared_vars->sims);
  NumericMatrix Pr0(shared_vars->n, shared_vars->sims);
  int cat_t = tt[0] ? shared_vars->cat_1 : shared_vars->cat_0;
  int cat_t_ctrl = tt[1] ? shared_vars->cat_1 : shared_vars->cat_0;
  int cat_c = tt[0] ? shared_vars->cat_0 : shared_vars->cat_1;
  int cat_c_ctrl = tt[1] ? shared_vars->cat_0 : shared_vars->cat_1;
  
  for(int j=0; j < shared_vars->sims; ++j){
    inner_loop(tt, j, Pr1, Pr0, cat_t, cat_t_ctrl, cat_c, cat_c_ctrl);
  }
  
  shared_vars->effects_tmp[e] = subtract_matrices(Pr1, Pr0);
  Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
}
//*/
void mediate_helper(Environment &env){
  Rcpp::Rcout<< "mediate_helper begin" << std::endl;
  shared_vars.reset(new shared_local_mediate_variables(env));
  /*
  for(std::size_t e=0;e<4;++e){
    outer_loop(e);
  }
   //*/
   /*
  env["et1"] = shared_vars->effects_tmp[0].to_Matrix();
  env["et2"] = shared_vars->effects_tmp[1].to_Matrix();
  env["et3"] = shared_vars->effects_tmp[2].to_Matrix();
  env["et4"] = shared_vars->effects_tmp[3].to_Matrix();
  //*/
  Rcpp::Rcout<< "mediate_helper end" << std::endl;
  shared_vars.reset();
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


shared_local_mediate_variables::shared_local_mediate_variables(Rcpp::Environment &env) : 
  n(env["n"]),
  sims(env["sims"]),
  cat_0(env["cat.0"]),
  cat_1(env["cat.1"]),
  treat(Rcpp::as<std::string>(env["treat"])),
  mediator(Rcpp::as<std::string>(env["mediator"])),
  PredictM0(DoubleMatrix::from_Rmatrix(Rcpp::as<NumericMatrix>(env["PredictM0"]))),
  PredictM1(DoubleMatrix::from_Rmatrix(Rcpp::as<NumericMatrix>(env["PredictM1"]))),
  YModel(DoubleMatrix::from_Rmatrix(Rcpp::as<NumericMatrix>(env["YModel"]))),
  y_data(DoubleMatrix::from_RDataFrame(Rcpp::as<DataFrame>(env["y.data"]))),
  tt_switch({{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}}})
  {
  effects_tmp.reserve(4);
  for(std::size_t i=0;i<4;++i){
    shared_vars->effects_tmp.emplace_back(n, sims);
  }
}
