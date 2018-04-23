// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>
#include <vector>
#include <nloptrAPI.h>

#include "cor_phylo.h"

using namespace Rcpp;


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




//' Inline C++ that does most of the work related to the log-likelihood function.
//' 
//' See below for the `nlopt` and `stats::optim` versions
//' 
//' 
//' @name cor_phylo_LL_
//' @noRd
//' 
//' 
inline double cor_phylo_LL_(const unsigned& n,
                            const arma::vec& par,
                            const arma::mat& XX,
                            const arma::mat& UU,
                            const arma::mat& MM,
                            const arma::mat& Vphy,
                            const arma::mat& tau,
                            const bool& REML,
                            const bool& constrain_d,
                            const bool& verbose) {
  
  
  uint_t n_ = Vphy.n_rows;
  uint_t p = XX.n_rows / n_;
  
  arma::mat L = make_L(par, n_, p);
  
  arma::mat R = L.t() * L;
  
  arma::vec d = make_d(par, n_, p, constrain_d, true);
  if (d.n_elem == 0) return MAX_RETURN;
  
  // OU transform
  arma::mat C = make_C(n_, p, tau, d, Vphy, R);
  
  arma::mat V = make_V(C, MM);
  double rcond_dbl = arma::rcond(V);
  if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
  
  arma::mat iV = arma::inv(V);
  arma::mat denom = tp(UU) * iV * UU;
  rcond_dbl = arma::rcond(denom);
  if (!arma::is_finite(rcond_dbl) || rcond_dbl < COND_MIN) return MAX_RETURN;
  
  arma::mat num = tp(UU) * iV * XX;
  arma::vec B0 = arma::solve(denom, num);
  arma::mat H = XX - UU * B0;
  
  double logdetV, det_sign;
  arma::log_det(logdetV, det_sign, iV);
  if (!arma::is_finite(logdetV)) return MAX_RETURN;
  logdetV *= -1;
  
  double LL;
  if (REML) {
    arma::mat to_det = tp(UU) * iV * UU;
    double det_val;
    arma::log_det(det_val, det_sign, to_det);
    double lhs = arma::as_scalar(tp(H) * iV * H);
    LL = 0.5 * (logdetV + det_val + lhs);
  } else {
    LL = 0.5 * arma::as_scalar(logdetV + tp(H) * iV * H);
  }
  
  if (verbose) {
    Rcout << LL << ' ';
    for (uint_t i = 0; i < par.n_elem; i++) Rcout << par(i) << ' ';
    Rcout << std::endl;
  }
  
  return LL;
}


//' `cor_phylo` log likelihood function for use with `nlopt`.
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
//'   In this function's case, it is an object of class `LL_info`.
//' 
//' @return the negative log likelihood
//' 
//' @name cor_phylo_LL
//' @noRd
//' 
double cor_phylo_LL(unsigned n, const double* x, double* grad, void* f_data) {
  
  LL_info* ll_info = (LL_info*) f_data;
  
  arma::vec par(n);
  for (uint i = 0; i < n; i++) par[i] = x[i];
  
  double LL = cor_phylo_LL_(n, par, ll_info->XX, ll_info->UU, ll_info->MM, 
                            ll_info->Vphy, ll_info->tau, ll_info->REML, 
                            ll_info->constrain_d, ll_info->verbose);
  
  ll_info->iters++;
  
  return LL;
}



//' `cor_phylo` log likelihood function for R's `stats::optim`.
//' 
//' 
//' @param par Initial values for the parameters to be optimized over.
//' @param ll_info_xptr `Rcpp::Xptr` object that points to a C++ `LL_info` object.
//'     This object stores all the other information needed for the log likelihood
//'     function.
//' 
//' @noRd
//' 
//' @name cor_phylo_LL_R
//' 
//[[Rcpp::export]]
double cor_phylo_LL_R(const arma::vec& par,
                      SEXP ll_info_xptr) {
  
  XPtr<LL_info> ll_info(ll_info_xptr);
  const unsigned n = par.n_elem;
  double LL = cor_phylo_LL_(n, par, ll_info->XX, ll_info->UU, ll_info->MM, 
                            ll_info->Vphy, ll_info->tau, ll_info->REML, 
                            ll_info->constrain_d, ll_info->verbose);
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
//' @inheritParams ll_info cp_get_output
//' @inheritParams max_iter cor_phylo
//' @inheritParams method cor_phylo
//' 
//' @return Nothing. `ll_info` is modified in place to have info from the model fit
//'   after this function is run.
//'
//' @name fit_cor_phylo_nlopt
//' @noRd
//' 
void fit_cor_phylo_nlopt(XPtr<LL_info> ll_info_xptr,
                         const double& rel_tol,
                         const int& max_iter,
                         const uint& method) {
  
  LL_info& ll_info(*ll_info_xptr);
  void* llop(&ll_info);
  
  unsigned n_pars = ll_info.par0.n_elem;
  
  // Dynamic array from ll_info
  // (from http://www.fredosaurus.com/notes-cpp/newdelete/50dynamalloc.html)
  double* x = NULL;   // initialize to nothing.
  x = new double[n_pars];  // Allocate n_pars doubles and save ptr in x.
  for (unsigned i = 0; i < n_pars; i++) x[i] = ll_info.par0(i);
  
  
  nlopt_algorithm alg;
  
  switch(method) {
  case 0: 
    alg = NLOPT_LN_NELDERMEAD;
    break;
  case 1: 
    alg = NLOPT_LN_BOBYQA;
    break;
  case 2: 
    alg = NLOPT_LN_SBPLX;
    break;
  case 3: 
    alg = NLOPT_LN_COBYLA;
    break;
  case 4: 
    alg = NLOPT_LN_PRAXIS;
    break;
  default:
    stop("Unknown method integer passed to fit_cor_phylo_nlopt");
  }
  
  nlopt_opt opt = nlopt_create(alg, n_pars);
  double min_LL;
  
  nlopt_set_min_objective(opt, cor_phylo_LL, llop);
  
  nlopt_set_ftol_rel(opt, rel_tol);
  nlopt_set_maxeval(opt, max_iter);
  
  nlopt_result res = nlopt_optimize(opt, x, &min_LL);
  ll_info.convcode = res;
  
  ll_info.LL = min_LL;
  for (unsigned i = 0; i < n_pars; i++) ll_info.min_par(i) = x[i];
  
  
  delete [] x;  // When done, free memory pointed to by x.
  x = NULL;     // Clear x to prevent using invalid memory reference.
  
  return;
}



//' Fit `cor_phylo` model using R's `stats::optim`.
//' 
//' Make sure this doesn't get run in parallel!
//'
//'
//' @inheritParams ll_info_xptr cor_phylo_LL_R
//' @inheritParams max_iter cor_phylo
//' @inheritParams method cor_phylo
//' 
//' @return Nothing. `ll_info_xptr` is modified in place to have info from the model fit
//'   after this function is run.
//'
//' @name fit_cor_phylo_R
//' @noRd
//' 
void fit_cor_phylo_R(XPtr<LL_info>& ll_info_xptr,
                     const double& rel_tol,
                     const int& max_iter,
                     const uint& method,
                     const std::vector<double>& sann) {
  
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  
  Rcpp::List opt;
  
  if (method > 5) {
    opt = optim(_["par"] = ll_info_xptr->par0,
                _["fn"] = Rcpp::InternalFunction(&cor_phylo_LL_R),
                _["ll_info"] = ll_info_xptr,
                _["method"] = "SANN",
                _["control"] = List::create(_["maxit"] = sann[0],
                                     _["temp"] = sann[1],
                                     _["tmax"] = sann[2],
                                     _["reltol"] = rel_tol));
    ll_info_xptr->par0 = as<arma::vec>(opt["par"]);
  }
  
  opt = optim(_["par"] = ll_info_xptr->par0,
              _["fn"] = Rcpp::InternalFunction(&cor_phylo_LL_R),
              _["ll_info"] = ll_info_xptr,
              _["method"] = "Nelder-Mead",
              _["control"] = List::create(_["maxit"] = max_iter,
                                   _["reltol"] = rel_tol));
  
  ll_info_xptr->min_par = as<arma::vec>(opt["par"]);
  
  ll_info_xptr->LL = as<double>(opt["value"]);
  ll_info_xptr->convcode = as<int>(opt["convergence"]);
  ll_info_xptr->iters = as<arma::vec>(opt["counts"])(0);
  
  if (ll_info_xptr->verbose) {
    Rcout << ll_info_xptr->LL << ' ';
    arma::vec& par(ll_info_xptr->min_par);
    for (uint_t i = 0; i < par.n_elem; i++) Rcout << par(i) << ' ';
    Rcout << std::endl;
  }
  
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




//' Make an `LL_info` object based on input matrices.
//' 
//' This `LL_info` is used for model fitting.
//' 
//' @inheritParams X cor_phylo_
//' @inheritParams U cor_phylo_
//' @inheritParams M cor_phylo_
//' @inheritParams Vphy_ cor_phylo_
//' @inheritParams REML_ cor_phylo_
//' @inheritParams constrain_d_ cor_phylo_
//' @inheritParams verbose_ cor_phylo_
//' 
//' @return a LL_info that contains info necessary for model fitting
//' 
//' @name LL_info
//' @noRd
//' 
LL_info::LL_info(const arma::mat& X,
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
  
  tau = arma::vec(n, arma::fill::ones) * Vphy.diag().t() - Vphy;
  
  arma::mat Xs = X;
  std::vector<arma::mat> Us = U;
  arma::mat Ms = M;
  standardize_matrices(Xs, Us, Ms);
  
  
  XX = arma::reshape(Xs, Xs.n_elem, 1);
  MM = flex_pow(Ms, 2);
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
        const arma::mat& x(Us[i]);
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
  
  par0 = arma::vec((static_cast<double>(p) / 2) * (1 + p) + p);
  par0.fill(0.5);
  for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
    par0(arma::span(j, k)) = L(arma::span(i, p-1), i);
    j = k + 1;
    k += (p - i - 1);
  }
  
  min_par = par0;
  
}



inline void main_output(arma::mat& corrs, arma::mat& B, arma::mat& B_cov, arma::vec& d,
                        const LL_info& ll_info, const uint_t& n, const uint_t& p,
                        const arma::mat& mean_sd_X, const std::vector<arma::vec>& sd_U) {
  
  arma::mat L = make_L(ll_info.min_par, n, p);
  
  arma::mat R = L.t() * L;
  
  corrs = make_corrs(R);
  
  d = make_d(ll_info.min_par, n, p, ll_info.constrain_d);
  
  // OU transform
  arma::mat C = make_C(n, p, ll_info.tau, d, ll_info.Vphy, R);
  
  arma::mat V = make_V(C, ll_info.MM);
  
  arma::mat iV = arma::inv(V);
  
  arma::mat denom = tp(ll_info.UU) * iV * ll_info.UU;
  
  arma::mat num = tp(ll_info.UU) * iV * ll_info.XX;
  arma::vec B0 = arma::solve(denom, num);
  
  make_B_B_cov(B, B_cov, B0, p, iV, ll_info.UU, mean_sd_X, sd_U);
  
  return;
}


//' Retrieve objects for output `cor_phylo` object.
//' 
//' @inheritParams X cor_phylo_
//' @inheritParams U cor_phylo_
//' @param ll_info an LL_info object that contains info necessary to fit the model.
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
                   XPtr<LL_info> ll_info_xptr,
                   const double& rel_tol,
                   const int& max_iter,
                   const uint& method,
                   const uint& boot,
                   const std::vector<double>& sann) {
  
  LL_info& ll_info(*ll_info_xptr);
  
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
  
  
  uint_t n = ll_info.Vphy.n_rows;
  uint_t p = mean_sd_X.n_rows;
  
  /*
   Get the main output from cor_phylo: correlations, coefficient estimates, 
   coefficient covariances, and phylogenetic signals.
   The objects constructed below correspond to this information.
   */
  arma::mat corrs;
  arma::mat B;
  arma::mat B_cov;
  arma::vec d;
  main_output(corrs, B, B_cov, d, ll_info, n, p, mean_sd_X, sd_U);
  
  double logLik = -0.5 * std::log(2 * arma::datum::pi);
  if (ll_info.REML) {
    logLik *= (n * p - ll_info.UU.n_cols);
    arma::mat to_det = tp(ll_info.XX) * ll_info.XX;
    double det_val, det_sign;
    arma::log_det(det_val, det_sign, to_det);
    logLik += 0.5 * det_val - ll_info.LL;
  } else {
    logLik *= (n * p);
    logLik -= ll_info.LL;
  }
  
  double k = ll_info.min_par.n_elem + ll_info.UU.n_cols;
  double AIC, BIC;
  AIC = -2 * logLik + 2 * k;
  BIC = -2 * logLik + k * std::log(n / arma::datum::pi);
  
  
  List boot_list = List::create();
  if (boot > 0) {
    // `boot_mats` stores matrices that we'll need for bootstrapping
    boot_mats bm(mean_sd_X, sd_U, B, d, ll_info, U);
    boot_results br(p, B.n_rows, boot);
    if (method < 5) {
      for (uint_t b = 0; b < boot; b++) {
        Rcpp::checkUserInterrupt();
        bm.one_boot(ll_info_xptr, br, b, rel_tol, max_iter, method);
      }
    } else {
      for (uint_t b = 0; b < boot; b++) {
        Rcpp::checkUserInterrupt();
        bm.one_boot(ll_info_xptr, br, b, rel_tol, max_iter, method, sann);
      }
    }
    std::vector<NumericMatrix> failed_mats(br.failed_inds.size());
    for (uint i = 0; i < br.failed_inds.size(); i++) {
      failed_mats[i] = wrap(br.failed_mats[i]);
    }
    boot_list = List::create(_["corrs"] = br.corrs, _["d"] = br.d, _["B0"] = br.B0,
                             _["B_cov"] = br.B_cov, _["failed"] = br.failed_inds,
                             _["failed_mats"] = failed_mats);
  }
  
  // Now the final output list
  List out = List::create(
    _["corrs"] = corrs,
    _["d"] = d,
    _["B"] = B,
    _["B_cov"] = B_cov,
    _["logLik"] = logLik,
    _["AIC"] = AIC,
    _["BIC"] = BIC,
    _["niter"] = ll_info.iters,
    _["convcode"] = ll_info.convcode,
    _["bootstrap"] = boot_list
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
//' @param method the `method` input to `cor_phylo`, converted to an integer
//'   that indexes which method it should be.
//' 
//' @return a list containing output information, to later be coerced to a `cor_phylo`
//'   object by the `cor_phylo` function.
//' @noRd
//' @name cor_phylo_
//' 
//[[Rcpp::export]]
List cor_phylo_(const arma::mat& X,
                const std::vector<arma::mat>& U,
                const arma::mat& M,
                const arma::mat& Vphy_,
                const bool& REML,
                const bool& constrain_d,
                const bool& verbose,
                const double& rel_tol,
                const int& max_iter,
                const uint& method,
                const uint& boot,
                const std::vector<double>& sann) {
  
  // LL_info is C++ class to use for organizing info for optimizing
  XPtr<LL_info> ll_info_xptr(new LL_info(X, U, M, Vphy_, REML, constrain_d, 
                                        verbose), true);

  // Do the fitting. `method < 5` means to use nlopt. Otherwise, use R's `stats::optim`.
  if (method < 5) {
    fit_cor_phylo_nlopt(ll_info_xptr, rel_tol, max_iter, method);
  } else {
    fit_cor_phylo_R(ll_info_xptr, rel_tol, max_iter, method, sann);
  }
  
  // Retrieve output from `ll_info` object and convert to list
  List output = cp_get_output(X, U, ll_info_xptr, rel_tol, max_iter, method, boot,
                              sann);
  
  return output;
  
}










/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Bootstrapping functions
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */



boot_mats::boot_mats(const arma::mat& mean_sd_X_, const std::vector<arma::vec>& sd_U_,
                     const arma::mat& B_, const arma::vec& d_, const LL_info& ll_info,
                     const std::vector<arma::mat>& U)
  : mean_sd_X(mean_sd_X_), sd_U(sd_U_), mean_sd_X0(mean_sd_X_),
    XX(ll_info.XX), MM(ll_info.MM), B0(B_.col(0)), iD(), U_add(), Us(U) {
  
  uint_t n = ll_info.Vphy.n_rows;
  uint_t p = mean_sd_X_.n_rows;
  
  arma::mat L = make_L(ll_info.min_par, n, p);
  arma::mat R = L.t() * L;
  arma::mat C = make_C(n, p, ll_info.tau, d_, ll_info.Vphy, R);
  arma::mat V = make_V(C, MM);
  iD = arma::chol(V).t();
  
  // For adding to X values
  if (ll_info.UU.n_cols > p) {
    U_add = ll_info.UU;
    U_add.cols(0, p).fill(0); // keeps means at zero;
    U_add = tp(B0) * tp(U_add);
    U_add = tp(U_add);
  } else {
    U_add = arma::mat(n * p, 1, arma::fill::zeros);
  }
  
  // Standardized U vector of matrices
  if (Us.size() > 0) {
    for (uint_t i = 0; i < U.size(); i++) {
      arma::mat& Ui(Us[i]);
      for (uint_t j = 0; j < Ui.n_cols; j++) {
        double sd = arma::stddev(Ui.col(j));
        Ui.col(j) -= arma::mean(Ui.col(j));
        if (sd > 0) Ui.col(j) /= sd;
      }
    }
  }
}

//' Iterate from a boot_mats object in prep for a bootstrap replicate.
//' 
//' This ultimately updates the LL_info object with new XX and MM matrices,
//' and updates the boot_results object with the mean and sd.
//' 
//' @param ll_info An LL_info object that will inform the next call to the
//'     log-likelihood function.
//' @param br A boot_results object that stores output from bootstrapping.
//' 
//' @name boot_mats_iterate
//' @noRd
//' 
void boot_mats::iterate(LL_info& ll_info) {
  
  // Reset this matrix to original
  mean_sd_X = mean_sd_X0;
  
  uint_t n = ll_info.Vphy.n_rows;
  uint_t p = mean_sd_X.n_cols;
  
  arma::mat M = MM;
  arma::vec rnd(rnorm(n * p));
  arma::mat X = iD * rnd + U_add;
  M.reshape(n, p);
  X.reshape(n, p);
  for (uint_t i = 0; i < p; i++) {
    double sd_ = arma::stddev(X.col(i));
    double mean_ = arma::mean(X.col(i));
    // Filling in mean_sd_X before normalizing
    // (I'm multiplying by original sd to convert to back non-normalized scale)
    mean_sd_X(i, 0) += (mean_ * mean_sd_X0(i, 1));
    mean_sd_X(i, 1) *= sd_;
    X.col(i) -= mean_;
    X.col(i) /= sd_;
    M.col(i) /= sd_;
  }
  
  arma::mat eps = X;
  
  arma::mat L;
  if (Us.size() > 0) {
    for (uint_t i = 0; i < p; i++) {
      if (Us[i].n_cols > 0) {
        const arma::mat& x(Us[i]);
        const arma::vec& y(X.col(i));
        arma::vec coef = arma::solve(x, y);
        arma::vec res = y - x * coef;
        eps.col(i) = res;
      } else {
        eps.col(i) = X.col(i) - arma::mean(X.col(i));
      }
    }
  }
  L = arma::cov(eps);
  L = arma::chol(L);
  L = tp(L);
  
  
  arma::vec& par0(ll_info.par0);
  
  par0 = arma::vec((static_cast<double>(p) / 2) * (1 + p) + p);
  par0.fill(0.5);
  for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
    par0(arma::span(j, k)) = L(arma::span(i, p-1), i);
    j = k + 1;
    k += (p - i - 1);
  }
  
  ll_info.min_par = par0;
  
  
  M.reshape(M.n_elem, 1);
  X.reshape(X.n_elem, 1);
  
  ll_info.XX = X;
  ll_info.MM = M;
  
  return;
}



// Method to return data when convergence fails
void boot_mats::failed(LL_info& ll_info, boot_results& br, const uint& i) {
  
  br.failed_inds.push_back(i+1);
  
  uint n = ll_info.Vphy.n_rows;
  uint p = mean_sd_X.n_rows;

  arma::mat X = ll_info.XX;
  arma::mat M = ll_info.MM;
  X.reshape(n, p);
  M.reshape(n, p);
  for (uint i = 0; i < p; i++) {
    X.col(i) *= mean_sd_X(i, 1);
    X.col(i) += mean_sd_X(i, 0);
    M.col(i) *= mean_sd_X(i, 1);
  }
  arma::mat XM = arma::join_horiz(X, M);
  br.failed_mats.push_back(XM);
  return;
}

void boot_mats::one_boot(XPtr<LL_info>& ll_info_xptr, boot_results& br,
                         const uint_t& i, const double& rel_tol, const int& max_iter,
                         const uint_t& method) {
  
  LL_info& ll_info(*ll_info_xptr);
  
  // Generate new data
  iterate(ll_info);
  
  // Do the fitting
  fit_cor_phylo_nlopt(ll_info_xptr, rel_tol, max_iter, method);
  
  if (ll_info_xptr->convcode < 0) failed(ll_info, br, i);
  
  arma::mat corrs;
  arma::mat B;
  arma::mat B_cov;
  arma::vec d;
  main_output(corrs, B, B_cov, d, ll_info, ll_info.Vphy.n_rows, mean_sd_X.n_rows,
              mean_sd_X, sd_U);
  
  // Add values to boot_results
  br.insert_values(i, corrs, d, B.col(0), B_cov);
  
  return;
}


void boot_mats::one_boot(XPtr<LL_info>& ll_info_xptr, boot_results& br,
                         const uint_t& i, const double& rel_tol, const int& max_iter,
                         const uint_t& method, const std::vector<double>& sann) {
  
  LL_info& ll_info(*ll_info_xptr);
  
  // Generate new data
  iterate(ll_info);
  
  // Do the fitting.
  fit_cor_phylo_R(ll_info_xptr, rel_tol, max_iter, method, sann);
  
  if (ll_info_xptr->convcode != 0) failed(ll_info, br, i);

  arma::mat corrs;
  arma::mat B;
  arma::mat B_cov;
  arma::vec d;
  main_output(corrs, B, B_cov, d, ll_info, ll_info.Vphy.n_rows, mean_sd_X.n_rows,
              mean_sd_X, sd_U);

  // Add values to boot_results
  br.insert_values(i, corrs, d, B.col(0), B_cov);
  
  return;
}

