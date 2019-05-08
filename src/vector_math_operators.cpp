#include "vector_math_operators.hpp"

// [[Rcpp::plugins(cpp11)]]

std::vector<double> wrap_DoubleMatrix_column_or_row(std::vector<double*> vec){
  std::vector<double> result;
  result.reserve(vec.size());
  for(std::size_t i=0;i<vec.size();++i){
    result.push_back(*vec[i]);
  }
  return result;
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