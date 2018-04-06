// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>
#include <vector>
#include <string>
#include <nloptrAPI.h>

#include "cor_phylo.h"

using namespace Rcpp;


#define VERY_SMALL_TOL 0.00000001490116119385
#define COND_MIN 0.0000000001
#define MAX_RETURN 10000000000






/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Log likelihood function
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */

//' Log likelihood function.
//' 
//' Note that this function is referred to the "objective function" in the `nlopt`
//' documentation and the input arguments should not be changed.
//' See
//' [here](https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#objective-function)
//' for more information.
//' 
//' @param n the number of optimization parameters
//' @param x an array of length `n` of the optimization parameters
//' @param grad an array that is not used here, but the `nlopt` documentation describes
//'   it as such: "an array of length `n` which should (upon return) be set to the 
//'   gradient of the function with respect to the optimization parameters at `x`."
//' @param f_data pointer to an object with additional information for the function.
//'   In this function's case, it is an object of class `LL_obj`.
//' 
//' @return the negative log likelihood
//' 
//' @name cor_phylo_LL
//' @noRd
//' 
double cor_phylo_LL(unsigned n, const double* x, double* grad, void* f_data) {
  
  LL_obj* ll_obj = (LL_obj*) f_data;
  
  arma::vec par(n);
  for (uint i = 0; i < n; i++) par[i] = x[i];
  
  uint_t n_ = ll_obj->Vphy.n_rows;
  uint_t p = ll_obj->XX.n_rows / n_;
  
  arma::mat L = make_L(par, n_, p);
  
  arma::mat R = L.t() * L;
  
  arma::vec d = make_d(par, n_, p, ll_obj->constrain_d, true);
  if (d.n_elem == 0) return MAX_RETURN;
  
  // OU transform
  arma::mat C = make_C(n_, p, ll_obj->tau, d, ll_obj->Vphy, R);
  
  arma::mat V = make_V(C, ll_obj->MM);
  double rcond_dbl = arma::rcond(V);
  if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
  
  arma::mat iV = arma::inv(V);
  arma::mat denom = tp(ll_obj->UU) * iV * ll_obj->UU;
  rcond_dbl = arma::rcond(denom);
  if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
  
  arma::mat num = tp(ll_obj->UU) * iV * ll_obj->XX;
  arma::vec B0 = arma::solve(denom, num);
  arma::mat H = ll_obj->XX - ll_obj->UU * B0;
  
  double logdetV, det_sign;
  arma::log_det(logdetV, det_sign, iV);
  if (!arma::is_finite(logdetV)) return MAX_RETURN;
  logdetV *= -1;
  
  double LL;
  if (ll_obj->REML) {
    arma::mat to_det = tp(ll_obj->UU) * iV * ll_obj->UU;
    double det_val;
    arma::log_det(det_val, det_sign, to_det);
    double lhs = arma::as_scalar(tp(H) * iV * H);
    LL = 0.5 * (logdetV + det_val + lhs);
  } else {
    LL = 0.5 * arma::as_scalar(logdetV + tp(H) * iV * H);
  }
  
  if (ll_obj->verbose) {
    Rcout << LL << ' ';
    for (uint_t i = 0; i < par.n_elem; i++) Rcout << par(i) << ' ';
    Rcout << std::endl;
  }
  
  return LL;
}






/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Fit using nlopt
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */


//' Fit cor_phylo model using nlopt.
//'
//'
//' @inheritParams ll_obj cp_get_output
//' @inheritParams max_iter cor_phylo
//' @inheritParams method cor_phylo
//' 
//' @return Nothing. `ll_obj` is modified in place to have info from the model fit
//'   after this function is run.
//'
//' @name fit_cor_phylo
//' @noRd
//' 
void fit_cor_phylo(LL_obj& ll_obj,
                   const int& max_iter,
                   const std::string& method) {
  
  void* llop(&ll_obj);
  
  unsigned n_pars = ll_obj.par0.n_elem;
  
  // Dynamic array from ll_obj
  // (from http://www.fredosaurus.com/notes-cpp/newdelete/50dynamalloc.html)
  double* x = NULL;   // initialize to nothing.
  x = new double[n_pars];  // Allocate n_pars doubles and save ptr in x.
  for (unsigned i = 0; i < n_pars; i++) x[i] = ll_obj.par0(i);
  
  nlopt_opt opt;

  if (method == "neldermead") {
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, n_pars);
  } else if (method == "sbplx") {
    opt = nlopt_create(NLOPT_LN_SBPLX, n_pars);
  } else if (method == "bobyqa") {
    opt = nlopt_create(NLOPT_LN_BOBYQA, n_pars);
  } else if (method == "cobyla") {
    opt = nlopt_create(NLOPT_LN_COBYLA, n_pars);
  } else if (method == "praxis") {
    opt = nlopt_create(NLOPT_LN_PRAXIS, n_pars);
  } else {
    stop("method not recognized. Use bobyqa, cobyla, or praxis.");
  }
  double min_LL;
  
  nlopt_set_min_objective(opt, cor_phylo_LL, llop);
  
  nlopt_set_ftol_rel(opt, VERY_SMALL_TOL);
  nlopt_set_maxeval(opt, max_iter);
  
  nlopt_result res = nlopt_optimize(opt, x, &min_LL);
  ll_obj.convcode = res;
  
  ll_obj.LL = min_LL;
  for (unsigned i = 0; i < n_pars; i++) ll_obj.min_par(i) = x[i];
  
  
  delete [] x;  // When done, free memory pointed to by x.
  x = NULL;     // Clear x to prevent using invalid memory reference.
  
  return;
}





/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Other functions
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */


//' Standardize matrices in place.
//' 
//' Makes each column of the `X` matrix have mean of zero and standard deviation of 1.
//' If `U` isn't empty, this function makes each column in each matrix have
//' mean of zero and standard deviation of 1, unless all values are the same, in which
//' case it keeps the standard deviation at zero.
//' Divides each column of `M` by the original standard deviation of that column in 
//' `X`.
//' 
//' 
//' @inheritParams X cor_phylo_
//' @inheritParams U cor_phylo_
//' @inheritParams M cor_phylo_
//' 
//' @return Nothing. Matrices are standardized in place.
//' 
//' @name standardize_matrices
//' @noRd
//' 
void standardize_matrices(arma::mat& X,
                          std::vector<arma::mat>& U,
                          arma::mat& M) {
  
  uint_t p = X.n_cols;
  
  for (uint_t i = 0; i < p; i++) {
    double sd = arma::stddev(X.col(i));
    X.col(i) -= arma::mean(X.col(i));
    X.col(i) /= sd;
    M.col(i) /= sd;
  }
  
  if (U.size() > 0) {
    for (uint_t i = 0; i < U.size(); i++) {
      arma::mat& Ui(U[i]);
      for (uint_t j = 0; j < Ui.n_cols; j++) {
        double sd = arma::stddev(Ui.col(j));
        Ui.col(j) -= arma::mean(Ui.col(j));
        if (sd > 0) Ui.col(j) /= sd;
      }
    }
  }
  
  return;
}




//' Make an `LL_obj` object based on input matrices.
//' 
//' This `LL_obj` is used for model fitting.
//' 
//' @inheritParams X cor_phylo_
//' @inheritParams U cor_phylo_
//' @inheritParams M cor_phylo_
//' @inheritParams Vphy_ cor_phylo_
//' @inheritParams REML_ cor_phylo_
//' @inheritParams constrain_d_ cor_phylo_
//' @inheritParams verbose_ cor_phylo_
//' 
//' @return a LL_obj that contains info necessary for model fitting
//' 
//' @name LL_obj
//' @noRd
//' 
LL_obj::LL_obj(const arma::mat& X,
               const std::vector<arma::mat>& U,
               const arma::mat& M,
               const arma::mat& Vphy_,
               const bool& REML_,
               const bool& constrain_d_,
               const bool& verbose_) 
  : REML(REML_), constrain_d(constrain_d_), verbose(verbose_), iters(0) {
  
  uint_t n = Vphy_.n_rows;
  uint_t p = X.n_cols;
  
  Vphy = Vphy_;
  Vphy /= Vphy_.max();
  double val, sign;
  arma::log_det(val, sign, Vphy);
  val = std::exp(val / n);
  Vphy /= val;
  
  arma::mat Xs = X;
  std::vector<arma::mat> Us = U;
  arma::mat SeMs = M;
  standardize_matrices(Xs, Us, SeMs);
  
  
  XX = arma::reshape(Xs, Xs.n_elem, 1);
  MM = flex_pow(SeMs, 2);
  MM.reshape(MM.n_elem, 1);
  UU = arma::kron(arma::eye<arma::mat>(p,p),
                  arma::mat(n, 1, arma::fill::ones));
  
  if (U.size() > 0) {
    arma::vec zeros(p, arma::fill::zeros);
    for (uint_t i = 0; i < p; i++) {
      arma::vec dd = zeros;
      dd[i] = 1;
      arma::mat u = arma::kron(dd, Us[i]);
      for (uint_t j = 0; j < u.n_cols; j++) {
        if (arma::diff(u.col(j)).max() > 0) {
          UU.insert_cols(UU.n_cols,1);
          UU.col(UU.n_cols-1) = u.col(j);
        }
      }
    }
  }
  
  arma::mat L;
  arma::mat eps = Xs;
  if (U.size() > 0) {
    for (uint_t i = 0; i < p; i++) {
      if (U[i].n_cols > 0) {
        arma::mat x = Us[i];
        const arma::vec& y(Xs.col(i));
        arma::vec coef = arma::solve(x, y);
        arma::vec res = y - x * coef;
        eps.col(i) = res;
      } else {
        eps.col(i) = Xs.col(i) - arma::mean(Xs.col(i));
      }
    }
  }
  L = arma::cov(eps);
  L = arma::chol(L);
  L = tp(L);
  
  tau = arma::vec(n, arma::fill::ones) * Vphy.diag().t() - Vphy;
  
  par0 = arma::vec((p / 2) * (1 + p) + p);
  par0.fill(0.5);
  for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
    par0(arma::span(j, k)) = L(arma::span(i, p-1), i);
    j = k + 1;
    k += (p - i - 1);
  }
  
  min_par = par0;
  
}



//' Retrieve objects for output `cor_phylo` object.
//' 
//' @inheritParams X cor_phylo_
//' @inheritParams U cor_phylo_
//' @param ll_obj an LL_obj object that contains info necessary to fit the model.
//'   After optimization, it contains info from the model fit.
//' 
//' @return a list containing output information, to later be coerced to a `cor_phylo`
//'   object by the `cor_phylo` function.
//' 
//' @name cp_get_output
//' @noRd
//' 
List cp_get_output(const arma::mat& X,
                   const std::vector<arma::mat>& U,
                   LL_obj& ll_obj) {
  
  /*
   Extract the only info we any longer need from X and U
   */
  arma::mat mean_sd_X(X.n_cols, 2);
  mean_sd_X.col(0) = arma::conv_to<arma::vec>::from(arma::mean(X));
  mean_sd_X.col(1) = arma::conv_to<arma::vec>::from(arma::stddev(X));
  std::vector<arma::vec> sd_U(U.size());
  for (uint i = 0; i < U.size(); i++) {
    if (U[i].n_cols > 0) sd_U[i] = arma::conv_to<arma::vec>::from(arma::stddev(U[i]));
  }
  
  
  uint_t n = ll_obj.Vphy.n_rows;
  uint_t p = ll_obj.XX.n_rows / n;
  
  arma::mat L = make_L(ll_obj.min_par, n, p);
  
  arma::mat R = L.t() * L;

  arma::mat corrs = make_corrs(R);
  
  arma::vec d = make_d(ll_obj.min_par, n, p, ll_obj.constrain_d);
  
  // OU transform
  arma::mat C = make_C(n, p, ll_obj.tau, d, ll_obj.Vphy, R);
  
  arma::mat V = make_V(C, ll_obj.MM);
  
  arma::mat iV = arma::inv(V);
  
  arma::mat denom = tp(ll_obj.UU) * iV * ll_obj.UU;
  
  arma::mat num = tp(ll_obj.UU) * iV * ll_obj.XX;
  arma::vec B0 = arma::solve(denom, num);
  
  arma::mat B;
  arma::vec sd_vec;
  arma::mat B_cov;
  make_sd_B_mat_cov(B, sd_vec, B_cov, B0, p, iV, ll_obj.UU, mean_sd_X, sd_U);
  
  double logLik = -0.5 * std::log(2 * arma::datum::pi);
  if (ll_obj.REML) {
    logLik *= (n * p - ll_obj.UU.n_cols);
    arma::mat to_det = tp(ll_obj.XX) * ll_obj.XX;
    double det_val, det_sign;
    arma::log_det(det_val, det_sign, to_det);
    logLik += 0.5 * det_val - ll_obj.LL;
  } else {
    logLik *= (n * p);
    logLik -= ll_obj.LL;
  }
  
  double k = ll_obj.min_par.n_elem + ll_obj.UU.n_cols;
  double AIC, BIC;
  AIC = -2 * logLik + 2 * k;
  BIC = -2 * logLik + k * std::log(n / arma::datum::pi);
  
  // /*
  //  Making `Rcpp::XPtr` smart pointers to store C++ objects in R because these 
  //  matrices shouldn't be needed directly in R.
  //  I like to convert these `XPtr` objects to `SEXP` to make it very explicit that
  //  I want to pass them back to R.
  //  */
  // // `cp_matrices` stores matrices that we'll need for bootstrapping
  // XPtr<cp_matrices> cpm(new cp_matrices(mean_sd_X, sd_U, ll_obj.XX, ll_obj.UU,
  //                                       ll_obj.MM, ll_obj.Vphy, R, V, C, B));
  // SEXP cpm_(cpm);
  // 
  // // Turning the `LL_obj` into a pointer
  // XPtr<LL_obj> llop(&ll_obj);
  // SEXP llop_(llop);
  
  
  // Now the final output list
  List out = List::create(
    _["corrs"] = corrs,
    _["d"] = d,
    _["B"] = B,
    _["B_cov"] = B_cov,
    _["logLik"] = logLik,
    _["AIC"] = AIC,
    _["BIC"] = BIC,
    _["niter"] = ll_obj.iters,
    _["convcode"] = ll_obj.convcode
  );
  
  return out;
}



//' Inner function to create necessary matrices and do model fitting.
//' 
//' @param X a n x p matrix with p columns containing the values for the n taxa.
//' @param U a list of p matrices corresponding to the p columns of `X`, with each 
//'   matrix containing independent variables for the corresponding column of `X`.
//' @param M a n x p matrix with p columns containing standard errors of the trait 
//'   values in `X`. 
//' @param Vphy_ phylogenetic variance-covariance matrix from the input phylogeny.
//' @inheritParams REML cor_phylo
//' @inheritParams constrain_d cor_phylo
//' @inheritParams verbose cor_phylo
//' @inheritParams max_iter cor_phylo
//' @inheritParams method cor_phylo
//' 
//' @return a list containing output information, to later be coerced to a `cor_phylo`
//'   object by the `cor_phylo` function.
//' @noRd
//' 
//[[Rcpp::export]]
List cor_phylo_(const arma::mat& X,
                const std::vector<arma::mat>& U,
                const arma::mat& M,
                const arma::mat& Vphy_,
                const bool& REML,
                const bool& constrain_d,
                const bool& verbose,
                const int& max_iter,
                const std::string& method) {
  
  // LL_obj is C++ class to use for nlopt optimizing
  LL_obj ll_obj(X, U, M, Vphy_, REML, constrain_d, verbose);
  
  // Do the fitting
  fit_cor_phylo(ll_obj, max_iter, method);
  
  // Retrieve output from `ll_obj` object and convert to list
  List output = cp_get_output(X, U, ll_obj);
  
  return output;
}
