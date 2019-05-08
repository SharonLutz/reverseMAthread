#ifndef CLASS_DOUBLEMATRIX_HPP
#define CLASS_DOUBLEMATRIX_HPP
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

class DoubleMatrix;
#include "double_matrix_indexed_vector_ref.hpp"

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

#endif //CLASS_DOUBLEMATRIX_HPP