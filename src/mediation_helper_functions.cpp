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
// [[Rcpp::export]]
NumericMatrix cpp_mult(NumericMatrix m1, NumericMatrix m2);
// [[Rcpp::export]]
NumericMatrix cpp_tmult(NumericMatrix m1, NumericMatrix m2);

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
  
  void set_col_names(const std::vector<std::string> &col_names);
  void set_row_names(const std::vector<std::string> &row_names);
  
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
  
  
  void assign_column(std::vector<double> col, std::size_t col_i);
  void assign_row(std::vector<double> row, std::size_t row_i);
  
  void assign_column(std::vector<double> col, std::string& col_name);
  void assign_row(std::vector<double> row, std::string& row_name);
  
  
  void assign_column(double val, std::size_t col_i);
  void assign_row(double val, std::size_t row_i);
  
  void assign_column(double val, std::string& col_name);
  void assign_row(double val, std::string& row_name);
  
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
  this->transposed = other.transposed;
  this->num_rows = other.num_rows;
  this->num_cols = other.num_cols;
  this->data = other.data;
  if(other.col_names_assigned){
    this->set_col_names(other.col_names);
  }
  if(other.row_names_assigned){
    this->set_row_names(other.row_names);
  }
  /*
  this->row_names = other.row_names;
  this->col_names = other.col_names;
  this->row_name_to_index = other.row_name_to_index;
  this->col_name_to_index = other.col_name_to_index;
  this->col_names_assigned = other.col_names_assigned;
  this->row_names_assigned = other.row_names_assigned;
  //*/
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
      row[col_i] = &(data[row_i*(this->num_cols) + col_i]);
    }
  }
  for(std::size_t col_i=0; col_i < this->num_cols; ++col_i){
    auto &col = access_via_cols[col_i];
    col.resize(this->num_rows, nullptr);
    for(std::size_t row_i=0; row_i < this->num_rows; ++row_i){
      col[row_i] = &(data[row_i*(this->num_cols) + col_i]);
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
  std::size_t real_index = this->transposed ? (col_i * row_count) + row_i : (row_i * col_count) + col_i;
  //Rcout << "using real index : " << real_index << " computed from [" << row_i << "," << col_i << ']' <<std::endl;
  return data[real_index];
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
  //Rcout << col_names << std::endl;
  for(long long i=0; i < col_names.size(); ++i){
    std::string col_name;
    col_name=col_names[i];
    names.push_back(col_name);
  }
  this->set_col_names(names);
}

void DoubleMatrix::set_row_names(Rcpp::CharacterVector row_names){
  std::vector<std::string> names;
  names.reserve(row_names.size());
  //Rcout << row_names << std::endl;
  for(long long i=0; i < row_names.size(); ++i){
    std::string row_name;
    row_name = row_names[i];
    names.push_back(row_name);
  }
  this->set_row_names(names);
}

void DoubleMatrix::set_col_names(const std::vector<std::string> &new_col_names){
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
      this->row_names.emplace_back(new_col_names[i]);
      this->row_name_to_index[new_col_names[i]] = i;
    }
    this->row_names_assigned = true;
  } else {
    this->col_names.clear();
    this->col_names.reserve(col_count);
    this->col_name_to_index.clear();
    for(std::size_t i=0; i < new_col_names.size(); ++i){
      this->col_names.emplace_back(new_col_names[i]);
      this->col_name_to_index[new_col_names[i]] = i;
    }
    this->col_names_assigned = true;
  }
}

void DoubleMatrix::set_row_names(const std::vector<std::string> &new_row_names){
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
      this->col_names.emplace_back(new_row_names[i]);
      this->col_name_to_index[new_row_names[i]] = i;
    }
    this->col_names_assigned = true;
  } else {
    this->row_names.clear();
    this->row_names.reserve(row_count);
    this->row_name_to_index.clear();
    for(std::size_t i=0; i < new_row_names.size(); ++i){
      this->row_names.emplace_back(new_row_names[i]);
      this->row_name_to_index[new_row_names[i]] = i;
    }
    this->row_names_assigned = true;
  }
}

std::string DoubleMatrix::get_row_name(std::size_t row_i){
  if(row_i >= this->get_num_rows()){
    std::stringstream ss;
    ss << "index out of bounds, " << row_i << " >= " << this->get_num_rows();
    Rf_error(ss.str().c_str());
  }
  if((this->transposed && this->col_names_assigned) || ((!this->transposed) && this->row_names_assigned) ){
    if(this->transposed){
      return this->col_names[row_i];
    } else {
      return this->row_names[row_i];
    }
  }else {
    std::stringstream ss;
    ss << row_i;
    return ss.str();
  }
}

std::string DoubleMatrix::get_col_name(std::size_t col_i){
  if(col_i >= this->get_num_cols()){
    std::stringstream ss;
    ss << "index out of bounds, " << col_i << " >= " << this->get_num_cols();
    Rf_error(ss.str().c_str());
  }
  if((this->transposed && this->row_names_assigned) || ((!this->transposed) && this->col_names_assigned) ){
    if(this->transposed){
      return this->row_names[col_i];
    } else {
      return this->col_names[col_i];
    }
  }else {
    std::stringstream ss;
    ss << "V" << col_i;
    return ss.str();
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
  
  Rcout << "dim1["<<rows_1<<','<<cols_1<<"]; dim2["<<rows_2<<','<<cols_2<<']'<<std::endl;
  
  return result;
  
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

void DoubleMatrix::assign_column(std::vector<double*> col, std::size_t col_i){
  std::vector<double *> target = this->operator[](col_i);
  if(col.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< col.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = *col[i];
  }
}

void DoubleMatrix::assign_row(std::vector<double*> row, std::size_t row_i){
  std::vector<double *> target = this->operator()(row_i);
  if(row.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< row.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = *row[i];
  }
}

void DoubleMatrix::assign_column(std::vector<double*> col, std::string& col_name){
  std::vector<double *> target = this->operator[](col_name);
  if(col.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< col.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = *col[i];
  }
}

void DoubleMatrix::assign_row(std::vector<double*> row, std::string& row_name){
  std::vector<double *> target = this->operator()(row_name);
  if(row.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< row.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = *row[i];
  }
}

void DoubleMatrix::assign_column(std::vector<double> col, std::size_t col_i){
  std::vector<double *> target = this->operator[](col_i);
  if(col.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< col.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = col[i];
  }
}

void DoubleMatrix::assign_row(std::vector<double> row, std::size_t row_i){
  std::vector<double *> target = this->operator()(row_i);
  if(row.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< row.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = row[i];
  }
}

void DoubleMatrix::assign_column(std::vector<double> col, std::string& col_name){
  std::vector<double *> target = this->operator[](col_name);
  if(col.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< col.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = col[i];
  }
}

void DoubleMatrix::assign_row(std::vector<double> row, std::string& row_name){
  std::vector<double *> target = this->operator()(row_name);
  if(row.size() != target.size()){
    std::stringstream ss;
    ss << "dimension mismatch : " << target.size() << " vs "<< row.size();
    Rf_error(ss.str().c_str());
  }
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = row[i];
  }
}


void DoubleMatrix::assign_column(double val, std::string& col_name){
  std::vector<double *> target = this->operator[](col_name);
  
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = val;
  }
}

void DoubleMatrix::assign_row(double val, std::string& row_name){
  std::vector<double *> target = this->operator()(row_name);
  
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = val;
  }
}

void DoubleMatrix::assign_column(double val, std::size_t col_i){
  std::vector<double *> target = this->operator[](col_i);
  
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = val;
  }
}

void DoubleMatrix::assign_row(double val, std::size_t row_i){
  std::vector<double *> target = this->operator()(row_i);
  
  for(std::size_t i=0;i<target.size();++i){
    *target[i] = val;
  }
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
    result(row_i, 0) = vec[row_i];
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_vector(std::vector<double> &svec){
  std::size_t num_rows = svec.size();
  DoubleMatrix result(num_rows, 1);
  for(unsigned long long row_i=0; row_i < num_rows;++row_i){
    result(row_i, 0) = svec[row_i];
  }
  return result;
}

DoubleMatrix DoubleMatrix::from_Rvector(Rcpp::NumericVector vec, std::size_t num_rows, std::size_t num_cols){
  std::size_t veclen = vec.length();
  if(num_rows * num_cols != veclen){
    std::stringstream ss;
    ss << "rows and columns do not match incoming vector length";
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
    ss << "rows and columns do not match incoming vector length";
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

template<int RTYPE>
int get_vector_type_impl(Vector<RTYPE> xin){
  return RTYPE;
}

int get_vector_type(SEXP xin){
  RCPP_RETURN_VECTOR(get_vector_type_impl, xin);
}

DoubleMatrix DoubleMatrix::from_RDataFrame(Rcpp::DataFrame df){
  std::size_t ncol = df.length();
  std::size_t nrow = df.nrows();
  DoubleMatrix result(nrow, ncol);
  CharacterVector col_names = df.attr("names");
  for(unsigned long long col_i=0; col_i < ncol;++col_i){
    Rcpp::NumericVector row_vec = Rcpp::as<Rcpp::NumericVector>(df[col_i]);
    for(unsigned long long row_i=0; row_i < nrow;++row_i){
      result(row_i,col_i) = row_vec[row_i];
    }
  }
  result.set_col_names(col_names);
  SEXP row_name_sexp = df.attr("row.names");
  CharacterVector char_vec;
  if(get_vector_type(row_name_sexp) == get_vector_type(char_vec)){
    CharacterVector row_names = row_name_sexp;
    result.set_row_names(row_names);
  }
  return result;
}

Rcpp::DataFrame DoubleMatrix::to_DataFrame(){
  Rcpp::DataFrame result;
  std::size_t n_cols = this->get_num_cols();
  std::size_t n_rows = this->get_num_rows();
  CharacterVector col_names(n_cols);
  for(std::size_t col_i =0; col_i < n_cols; ++col_i){
    Rcpp::NumericVector col_vec(n_rows);
    std::string col_name = this->get_col_name(col_i);
    //Rcout << "col_name[" << col_i << "]:" <<col_name<< std::endl;
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      col_vec[row_i] = this->operator()(row_i, col_i);
    }
    result[col_name]=col_vec;
  }
  if(this->row_names_assigned){
    //Rcout << "row names assigned" << std::endl;
    CharacterVector row_names(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      row_names[row_i] = this->get_row_name(row_i);
    }
    result.attr("row.names") = row_names;
  } else {
    //Rcout << "row names not assigned" << std::endl;
    IntegerVector row_names(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      row_names[row_i]=row_i;
    }
    result.attr("row.names") = row_names;
  }
  result.attr("class") = "data.frame";
  return result;
}

Rcpp::NumericMatrix DoubleMatrix::to_Matrix(){
  Rcpp::NumericMatrix result(this->get_num_rows(), this->get_num_cols());
  for(std::size_t row_i=0;row_i < this->get_num_rows(); ++row_i){
    for(std::size_t col_i=0;col_i < this->get_num_cols(); ++col_i){
      result(row_i,col_i) = this->operator()(row_i, col_i);
    }
  }
  if(this->row_names_assigned){
    if(this->transposed){
      colnames(result) = wrap(this->row_names);
    } else {
      rownames(result) = wrap(this->row_names);
    }
  }
  if(this->col_names_assigned){
    if(this->transposed){
      rownames(result) = wrap(this->col_names);
    } else {
      colnames(result) = wrap(this->col_names);
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
  std::vector<std::string> terms;
  std::array< std::array<int, 4>, 4 > tt_switch;
  DoubleMatrix PredictM0;
  DoubleMatrix PredictM1;
  DoubleMatrix YModel;
  DoubleMatrix y_data;
  std::vector<DoubleMatrix>effects_tmp;
  shared_local_mediate_variables();
  static shared_local_mediate_variables from_environment(Rcpp::Environment & env);
};

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
shared_local_mediate_variables::shared_local_mediate_variables() : 
  PredictM0(0,0),
  PredictM1(0,0),
  YModel(0,0),
  y_data(0,0){}

shared_local_mediate_variables shared_local_mediate_variables::from_environment(Rcpp::Environment & env) {
  shared_local_mediate_variables result;
  result.tt_switch = {{{1,1,1,0},{0,0,1,0},{1,0,1,1},{1,0,0,0}}};
  result.n = env["n"];
  result.sims = env["sims"];
  result.cat_0 = env["cat.0"];
  result.cat_1 = env["cat.1"];
  result.treat = Rcpp::as<std::string>(env["treat"]);
  result.mediator = Rcpp::as<std::string>(env["mediator"]);
  
  NumericMatrix R_PredictM0 = env["PredictM0"];
  NumericMatrix R_PredictM1 = env["PredictM1"];
  NumericMatrix R_YModel = env["YModel"];
  SEXP df_data = env["y.data"];
  DataFrame R_y_data = df_data;
  
  Rcpp::Language terms = R_y_data.attr("terms");
  CharacterVector term_labels = terms.attr("term.labels");
  
  for(long long int i=0;i<term_labels.size();++i){
    std::string label;
    label = term_labels[i];
    result.terms.emplace_back(label);
  }
  
  result.PredictM0 = DoubleMatrix::from_Rmatrix(R_PredictM0);
  result.PredictM1 = DoubleMatrix::from_Rmatrix(R_PredictM1);
  
  result.YModel = DoubleMatrix::from_Rmatrix(R_YModel);
  result.y_data = DoubleMatrix::from_RDataFrame(R_y_data);
  
  result.effects_tmp.reserve(4);
  for(std::size_t i=0;i<4;++i){
    result.effects_tmp.emplace_back(result.n, result.sims);
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

//operator*
std::vector<double> operator*(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result = v1;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] * v2[i];
  }
  return result;
}

std::vector<double> operator*(std::vector<double> v, double val){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] * val;
  }
  return result;
}

std::vector<double> operator*(double val, std::vector<double> v){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] * val;
  }
  return result;
}
//operator+
std::vector<double> operator+(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result = v1;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] + v2[i];
  }
  return result;
}

std::vector<double> operator+(std::vector<double> v, double val){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] + val;
  }
  return result;
}

std::vector<double> operator+(double val, std::vector<double> v){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] + val;
  }
  return result;
}
//operator-
std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result = v1;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] - v2[i];
  }
  return result;
}

std::vector<double> operator-(std::vector<double> v, double val){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] - val;
  }
  return result;
}

std::vector<double> operator-(double val, std::vector<double> v){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = val - result[i];
  }
  return result;
}
//operator/
std::vector<double> operator/(std::vector<double> v1, std::vector<double> v2){
  std::vector<double> result = v1;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] / v2[i];
  }
  return result;
}

std::vector<double> operator/(std::vector<double> v, double val){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = result[i] / val;
  }
  return result;
}

std::vector<double> operator/(double val, std::vector<double> v){
  std::vector<double> result=v;
  for(std::size_t i=0;i<result.size();++i){
    result[i] = val / result[i];
  }
  return result;
}

DoubleMatrix pred_to_model_mat(shared_local_mediate_variables& sv, DoubleMatrix &pred_mat){
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

std::vector<double> wrap_DoubleMatrix_column_or_row(std::vector<double*> vec){
  std::vector<double> result;
  result.reserve(vec.size());
  for(std::size_t i=0;i<vec.size();++i){
    result.push_back(*vec[i]);
  }
  return result;
}

//*
void inner_loop(
    shared_local_mediate_variables& sv,
    std::array<int, 4>& tt,
std::size_t j, 
DoubleMatrix &Pr1, 
DoubleMatrix &Pr0,
int cat_t,
int cat_t_ctrl,
int cat_c,
int cat_c_ctrl
) {
Rcpp::Rcout << "inner: " << j << std::endl;

DoubleMatrix pred_data_t = sv.y_data;

DoubleMatrix pred_data_c = sv.y_data;

Rcpp::Rcout << "TEST0" << std::endl;

pred_data_t.assign_column(cat_t, sv.treat);
pred_data_c.assign_column(cat_c, sv.treat);

Rcpp::Rcout << "TEST1" << std::endl;

//PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
//PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])

double tt_val_2 = tt[2];
double tt_val_3 = tt[3];

Rcpp::Rcout << "sv.PredictM1(j) size: " << sv.PredictM1(j).size() << std::endl;
Rcpp::Rcout << "sv.PredictM0(j) size: " << sv.PredictM0(j).size() << std::endl;

std::vector<double> PredictM1_row_j = wrap_DoubleMatrix_column_or_row(sv.PredictM1(j));
std::vector<double> PredictM0_row_j = wrap_DoubleMatrix_column_or_row(sv.PredictM0(j));

Rcpp::Rcout << "PredictM1_row_j size: " << PredictM1_row_j.size() << std::endl;
Rcpp::Rcout << "PredictM0_row_j size: " << PredictM0_row_j.size() << std::endl;

std::vector<double> PredictMt = (PredictM1_row_j * tt_val_2) + (PredictM0_row_j * (1-tt_val_2));
std::vector<double> PredictMc = (PredictM1_row_j * tt_val_3) + (PredictM0_row_j * (1-tt_val_3));
//NumericVector PredictMt = (sv.PredictM1(j,_) * tt[2]) + (sv.PredictM0(j,_) * (1 - tt[2]));
//NumericVector PredictMc = (sv.PredictM1(j,_) * tt[3]) + (sv.PredictM0(j,_) * (1 - tt[3]));

Rcpp::Rcout << "PredictMt size: " << PredictMt.size() << std::endl;
Rcpp::Rcout << "PredictMc size: " << PredictMc.size() << std::endl;

Rcpp::Rcout << "TEST2" << std::endl;

//pred.data.t[,mediator] <- PredictMt
//pred.data.c[,mediator] <- PredictMc
pred_data_t.assign_column(PredictMt, sv.mediator);
pred_data_c.assign_column(PredictMc, sv.mediator);

Rcpp::Rcout << "TEST3" << std::endl;

DoubleMatrix ymat_t = pred_to_model_mat(sv, pred_data_t);
DoubleMatrix ymat_c = pred_to_model_mat(sv, pred_data_c);

Rcpp::Rcout << "TEST4" << std::endl;

std::vector<double> YModel_row_j =  wrap_DoubleMatrix_column_or_row(sv.YModel(j));

DoubleMatrix YModel_j_as_matrix = DoubleMatrix::from_vector(YModel_row_j);

Rcpp::Rcout << "TEST5" << std::endl;

YModel_j_as_matrix.transpose();
ymat_t.transpose();
ymat_c.transpose();

Rcpp::Rcout << "TEST6" << std::endl;

DoubleMatrix mmult1 = YModel_j_as_matrix * ymat_t;
DoubleMatrix mmult2 = YModel_j_as_matrix * ymat_t;

Rcpp::Rcout << "TEST7" << std::endl;

Rcout << mmult1.get_num_rows() << ',' << mmult1.get_num_cols() << std::endl;
Rcout << mmult2.get_num_rows() << ',' << mmult2.get_num_cols() << std::endl;
//Pr1(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_t);
//Pr0(_,j) = cpp_matrix_mult_of_transposed(YModel_j_as_matrix,ymat_c);
}

void outer_loop(shared_local_mediate_variables& sv, std::size_t e) {
Rcpp::Rcout<< "outer_loop begin: " << e << std::endl;
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

sv.effects_tmp[e] = Pr1 - Pr0;
Rcpp::Rcout<< "outer_loop end: " << e << std::endl;
}
//*/


void mediate_helper(Environment &env){
  Rcpp::Rcout<< "mediate_helper begin" << std::endl;
  //shared_vars.reset(new shared_local_mediate_variables(env));
  shared_local_mediate_variables shared_vars = shared_local_mediate_variables::from_environment(env);
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
  Rcpp::Rcout<< "mediate_helper end" << std::endl;
  //shared_vars.reset();
}

NumericMatrix cpp_mult(NumericMatrix m1, NumericMatrix m2){
  DoubleMatrix dm1 = DoubleMatrix::from_Rmatrix(m1);
  DoubleMatrix dm2 = DoubleMatrix::from_Rmatrix(m2);
  return (dm1 * dm2).to_Matrix();
}

NumericMatrix cpp_tmult(NumericMatrix m1, NumericMatrix m2){
  DoubleMatrix dm1 = DoubleMatrix::from_Rmatrix(m1);
  DoubleMatrix dm2 = DoubleMatrix::from_Rmatrix(m2);
  dm1.transpose();
  dm2.transpose();
  return (dm1 * dm2).to_Matrix();
}
