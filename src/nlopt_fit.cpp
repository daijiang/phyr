// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


using namespace Rcpp;


#include <vector>
#include <nloptrAPI.h>

// This is the default relative tolerance from stats::optim() in R:
#define VERY_SMALL_TOL 0.00000001490116119385

using namespace Rcpp;


// Example structure to hold a matrix as input data and an unsigned integer to hold
// the number of number of times it took to find the minimum value
struct my_data {
  arma::mat xy;
  unsigned count;
  my_data(const arma::mat& m) 
    : xy(m), count(0) {};
};


// Example function to minimize over: sum of squares for a simple y = b0 + b1 * x
// regression
double sum_of_squares(unsigned n, const double* x, double* grad, void* f_data) {
  
  // Convert void pointer to the data structure I defined above:
  my_data* d = (my_data*)f_data;
  
  // Count sum of squares
  double out = 0;
  for (unsigned i = 0; i < d->xy.n_rows; i++) {
    double r = d->xy(i,1) - (x[0] + x[1] * d->xy(i,0));
    out += (r * r);
  }
  // Update count variable in the my_data object:
  d->count++;
  
  return out;
}





//' Testing optimization using nlopt: least squares for a simple regression
//' 
//' @param xy A two-column matrix with the X variable as the first column and the Y
//'     variable as the second one.
//' @param b0 First estimate of the beta0 regression estimate.
//' @param b1 First estimate of the beta1 regression estimate.
//' @param max_iter Maximum calls to the optimization function.
//' @param method Algorithm used for optimization. For now, the options are
//'     "bobyqa", "cobyla", or "praxis".
//'     See \url{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms}
//'     for more information on these algorithms.
//' 
//' @return A numeric vector with best set of parameters (b0 then b1) and
//'     sum of squares corresponding to that set.
//' 
//' @export
//' 
//[[Rcpp::export]]
std::vector<double> test_nlopt(const arma::mat& xy,
                               const double& b0,
                               const double& b1,
                               const int& max_iter,
                               const char* method) {
  
  if (max_iter < 0) stop("max_iter cannot be less than 0.");
  if (xy.n_cols != 2) stop("xy must have exactly two columns.");
  
  my_data md(xy);
  void* mdp(&md);
  
  unsigned n = 2;
  // initial guesses
  double x[2] = { b0, b1 };
  double min_f; // the minimum objective value, upon return
  
  nlopt_opt opt;
  if (strcmp(method, "bobyqa") == 0) {
    opt = nlopt_create(NLOPT_LN_BOBYQA, n);
  } else if (strcmp(method, "cobyla") == 0) {
    opt = nlopt_create(NLOPT_LN_COBYLA, n);
  } else if (strcmp(method, "praxis") == 0) {
    opt = nlopt_create(NLOPT_LN_PRAXIS, n);
  } else {
    stop("method not recognized. Use bobyqa, cobyla, or praxis.");
  }
  
  nlopt_set_min_objective(opt, sum_of_squares, mdp);
  
  nlopt_set_ftol_rel(opt, VERY_SMALL_TOL);
  nlopt_set_maxeval(opt, max_iter);
  
  nlopt_result res = nlopt_optimize(opt, x, &min_f);
  if (res < 0) {
    Rcpp::Rcout << "nlopt failed!\\n";
  } else {
    Rcout << "Number of iterations: " << md.count << std::endl;
  }
  
  std::vector<double> result(3);
  result[0] = x[0];
  result[1] = x[1];
  result[2] = min_f;
  return result;
}
