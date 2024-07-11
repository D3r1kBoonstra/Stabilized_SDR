#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate determinant of a 2x2 matrix
double determinant(const NumericMatrix& mat) {
  return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
}

// Function to calculate quadratic discriminant
NumericMatrix qda_shrink(NumericMatrix x, IntegerVector grouping) {
  // Convert x to matrix
  NumericMatrix data(x);
  
  // Get unique classes
  IntegerVector classes = unique(grouping);
  
  // Prepare storage
  int p = x.ncol();
  int N = x.nrow();
  int n_classes = classes.size();
  std::vector<int> n(n_classes);
  std::vector<double> priors(n_classes);
  std::vector<NumericVector> xbar(n_classes);
  std::vector<NumericMatrix> S(n_classes);
  std::vector<NumericMatrix> S_inv(n_classes);
  std::vector<NumericVector> out(N);
  
  // Calculate statistics for each class
  for (int i = 0; i < n_classes; ++i) {
    // Extract data for the current class
    NumericMatrix class_data = data[which(grouping == classes[i]), _];
    int class_size = class_data.nrow();
    n[i] = class_size;
    
    // Calculate prior probability
    priors[i] = static_cast<double>(class_size) / N;
    
    // Calculate mean
    xbar[i] = colMeans(class_data);
    
    // Calculate covariance matrix
    S[i] = cov(class_data);
    
    // Calculate inverse of covariance matrix
    S_inv[i] = solve(S[i]);
  }
  
  // Calculate discriminant function
  for (int k = 0; k < n_classes; ++k) {
    NumericMatrix d(N, 1);
    for (int i = 0; i < N; ++i) {
      NumericVector diff = as<NumericVector>(x.row(i)) - xbar[k];
      d(i, 0) = -0.5 * (sum(diff * (S_inv[k] * diff)) - determinant(S_inv[k])) + log(priors[k]);
    }
    out[k] = d;
  }
  
  // Find class with maximum discriminant for each observation
  NumericMatrix result(N, 1);
  for (int i = 0; i < N; ++i) {
    int max_index = 0;
    double max_value = out[0][i];
    for (int k = 1; k < n_classes; ++k) {
      if (out[k][i] > max_value) {
        max_index = k;
        max_value = out[k][i];
      }
    }
    result(i, 0) = classes[max_index];
  }
  
  return result;
}

/*** R
# To expose this function to R, you can use the following Rcpp function:
# sourceCpp("filename.cpp")
*/
