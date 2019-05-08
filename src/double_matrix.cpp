#include "double_matrix.hpp"

#include "vector_math_operators.hpp"

// [[Rcpp::plugins(cpp11)]]

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
  
  // Rcpp::Rcout << "dim1["<<rows_1<<','<<cols_1<<"]; dim2["<<rows_2<<','<<cols_2<<']'<<std::endl;
  
  for(unsigned long long int result_row_i=0;result_row_i<rows_1;++result_row_i){
    
    auto &row_1_i = this->operator()(result_row_i);
    
    for(unsigned long long int result_col_i=0; result_col_i<cols_2;++result_col_i){
      
      auto &col_2_i = other.operator[](result_col_i);
      
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
  bool rows_mismatch = this->get_num_rows() != other.get_num_rows();
  bool cols_mismatch = this->get_num_cols() != other.get_num_cols();
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
  bool rows_mismatch = this->get_num_rows() != other.get_num_rows();
  bool cols_mismatch = this->get_num_cols() != other.get_num_cols();
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
  Rcpp::IntegerVector dim = mat.attr("dim");
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
int get_vector_type_impl(Rcpp::Vector<RTYPE> xin){
  return RTYPE;
}

int get_vector_type(SEXP xin){
  RCPP_RETURN_VECTOR(get_vector_type_impl, xin);
}

DoubleMatrix DoubleMatrix::from_RDataFrame(Rcpp::DataFrame df){
  std::size_t ncol = df.length();
  std::size_t nrow = df.nrows();
  DoubleMatrix result(nrow, ncol);
  Rcpp::CharacterVector col_names = df.attr("names");
  for(unsigned long long col_i=0; col_i < ncol;++col_i){
    Rcpp::NumericVector row_vec = Rcpp::as<Rcpp::NumericVector>(df[col_i]);
    for(unsigned long long row_i=0; row_i < nrow;++row_i){
      result(row_i,col_i) = row_vec[row_i];
    }
  }
  result.set_col_names(col_names);
  SEXP row_name_sexp = df.attr("row.names");
  Rcpp::CharacterVector char_vec;
  if(get_vector_type(row_name_sexp) == get_vector_type(char_vec)){
    Rcpp::CharacterVector row_names = row_name_sexp;
    result.set_row_names(row_names);
  }
  return result;
}

Rcpp::DataFrame DoubleMatrix::to_DataFrame(){
  Rcpp::DataFrame result;
  std::size_t n_cols = this->get_num_cols();
  std::size_t n_rows = this->get_num_rows();
  Rcpp::CharacterVector col_names(n_cols);
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
    Rcpp::CharacterVector row_names(n_rows);
    for(std::size_t row_i =0; row_i < n_rows; ++row_i){
      row_names[row_i] = this->get_row_name(row_i);
    }
    result.attr("row.names") = row_names;
  } else {
    //Rcout << "row names not assigned" << std::endl;
    Rcpp::IntegerVector row_names(n_rows);
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
      colnames(result) = Rcpp::wrap(this->row_names);
    } else {
      rownames(result) = Rcpp::wrap(this->row_names);
    }
  }
  if(this->col_names_assigned){
    if(this->transposed){
      rownames(result) = Rcpp::wrap(this->col_names);
    } else {
      colnames(result) = Rcpp::wrap(this->col_names);
    }
  }
  return result;
}
