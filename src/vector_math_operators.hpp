#ifndef VECTOR_MATH_OPERATORS_HPP
#define VECTOR_MATH_OPERATORS_HPP
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

std::vector<double> wrap_DoubleMatrix_column_or_row(std::vector<double*> vec);

std::vector<double> operator*(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator*(std::vector<double> v, double val);
std::vector<double> operator*(double val, std::vector<double> v);

std::vector<double> operator+(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator+(std::vector<double> v, double val);
std::vector<double> operator+(double val, std::vector<double> v);

std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator-(std::vector<double> v, double val);
std::vector<double> operator-(double val, std::vector<double> v);

std::vector<double> operator/(std::vector<double> v1, std::vector<double> v2);
std::vector<double> operator/(std::vector<double> v, double val);
std::vector<double> operator/(double val, std::vector<double> v);

#endif //VECTOR_MATH_OPERATORS_HPP