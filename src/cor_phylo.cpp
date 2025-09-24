// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/*
 Prevents the warning "solve(): system is singular; attempting approx solution"
 */
#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>
#include <vector>

#include "cor_phylo.h"

using namespace Rcpp;


#define MAX_RETURN 10000000000.0






/*
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 
 Log likelihood function
 
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */


// `cor_phylo` log likelihood function.
// 
//[[Rcpp::export]]
double cor_phylo_LL(NumericVector par,
                    const arma::mat& XX,
                    const arma::mat& UU,
                    const arma::mat& MM,
                    const arma::mat& Vphy,
                    const arma::mat& tau,
                    const bool& REML,
                    const bool& constrain_d,
                    const double& lower_d,
                    const bool& verbose,
                    const double& rcond_threshold) {

  bool return_max = false;
  
  uint_t n = Vphy.n_rows;
  uint_t p = XX.n_rows / n;
  
  arma::mat L = make_L(par, p);
  
  arma::mat R = L.t() * L;
  
  arma::vec d = make_d(par, p, constrain_d, lower_d, return_max);
  if (return_max) return MAX_RETURN;
  
  // OU transform
  arma::mat C = make_C(n, p, tau, d, Vphy, R);
  
  arma::mat V = make_V(C, MM);
  if (V.has_nan() || V.has_inf()) return MAX_RETURN;
  double rcond_dbl = 0;
  rcond_dbl = arma::rcond(V);
  if (!std::isfinite(rcond_dbl) || rcond_dbl < rcond_threshold) return MAX_RETURN;
  
  arma::mat iV = arma::inv(V);
  arma::mat denom = UU.t() * iV * UU;
  if (denom.has_nan() || denom.has_inf()) return MAX_RETURN;
  rcond_dbl = arma::rcond(denom);
  if (!std::isfinite(rcond_dbl) || rcond_dbl < rcond_threshold) return MAX_RETURN;
  
  arma::mat num = UU.t() * iV * XX;
  arma::vec B0 = arma::solve(denom, num);
  arma::mat H = XX - UU * B0;
  
  double logdetV, det_sign;
  arma::log_det(logdetV, det_sign, iV);
  if (!std::isfinite(logdetV)) return MAX_RETURN;
  logdetV *= -1;
  
  double LL;
  if (REML) {
    arma::mat to_det = UU.t() * iV * UU;
    double det_val;
    arma::log_det(det_val, det_sign, to_det);
    double lhs = arma::as_scalar(H.t() * iV * H);
    LL = 0.5 * (logdetV + det_val + lhs);
  } else {
    LL = 0.5 * arma::as_scalar(logdetV + H.t() * iV * H);
  }
  
  if (verbose) {
    Rcout << LL << ' ';
    for (uint_t i = 0; i < static_cast<uint_t>(par.size()); i++) Rcout << par[i] << ' ';
    Rcout << std::endl;
  }
  
  return LL;
}







/*
 Return reciprocal condition numbers for matrices in the log likelihood function.
 
 This function is largely a repeat of the first part of the likelihood function.
 It is used in the output to guide users wanting to change the `rcond_threshold`
 argument.
 */
std::vector<double> return_rcond_vals(const LogLikInfo& ll_info) {
  
  const arma::vec& par(ll_info.min_par);
  const arma::mat& XX(ll_info.XX);
  const arma::mat& UU(ll_info.UU);
  const arma::mat& MM(ll_info.MM);
  const arma::mat& Vphy(ll_info.Vphy);
  const arma::mat& tau(ll_info.tau);
  // const bool& REML(ll_info.REML);
  const bool& constrain_d(ll_info.constrain_d);
  const double& lower_d(ll_info.lower_d);
  // const bool& verbose(ll_info.verbose);
  // const double& rcond_threshold(ll_info.rcond_threshold);
  
  std::vector<double> rconds_out(2);
  
  uint_t n = Vphy.n_rows;
  uint_t p = XX.n_rows / n;
  
  arma::mat L = make_L(par, p);
  
  arma::mat R = L.t() * L;
  
  arma::vec d = make_d(par, p, constrain_d, lower_d);

  // OU transform
  arma::mat C = make_C(n, p, tau, d, Vphy, R);
  
  arma::mat V = make_V(C, MM);
  double rcond_dbl = 0;
  if (V.has_nan() || V.has_inf()) {
    rcond_dbl = arma::datum::nan;
  } else {
    rcond_dbl = arma::rcond(V);
  }
  rconds_out[0] = rcond_dbl;
  
  arma::mat iV = arma::inv(V);
  arma::mat denom = UU.t() * iV * UU;
  if (denom.has_nan() || denom.has_inf()) {
    rcond_dbl = arma::datum::nan;
  } else {
    rcond_dbl = arma::rcond(denom);
  }
  rconds_out[1] = rcond_dbl;
  
  return rconds_out;
}







/*
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 
 Fit using nlopt
 
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */

/*
 Fit cor_phylo model using nlopt.
 */
void fit_cor_phylo_nlopt(LogLikInfo& ll_info,
                         const double& rel_tol,
                         const int& max_iter,
                         const std::string& method) {
  
  Rcpp::Environment nloptr_pkg = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptr_pkg["nloptr"];
  
  Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
  Rcpp::Function cor_phylo_LL_fxn = phyr_pkg["cor_phylo_LL"];
  
  std::string nlopt_algor;
  
  if (method == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD";
  if (method == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA";
  if (method == "subplex") nlopt_algor = "NLOPT_LN_SBPLX";
  
  List options = List::create(_["algorithm"] = nlopt_algor,
                              _["ftol_rel"] = rel_tol,
                              _["ftol_abs"] = rel_tol,
                              _["xtol_rel"] = 0.0001,
                              _["maxeval"] = max_iter);
  
  NumericVector par0(ll_info.par0.begin(), ll_info.par0.end());
  
  List opt = nloptr(_["x0"] = par0,
                   _["eval_f"] = cor_phylo_LL_fxn,
                   _["opts"] = options,
                   _["XX"] = ll_info.XX,
                   _["UU"] = ll_info.UU,
                   _["MM"] = ll_info.MM,
                   _["Vphy"] = ll_info.Vphy,
                   _["tau"] = ll_info.tau,
                   _["REML"] = ll_info.REML,
                   _["constrain_d"] = ll_info.constrain_d,
                   _["lower_d"] = ll_info.lower_d,
                   _["verbose"] = ll_info.verbose,
                   _["rcond_threshold"] = ll_info.rcond_threshold);
  
  ll_info.min_par = as<arma::vec>(opt["solution"]);
  
  ll_info.LL = as<double>(opt["objective"]);
  int convcode_ = as<int>(opt["status"]);
  
  if (convcode_ > 0) {
    if (convcode_ < 5) {
      ll_info.convcode = 0;
    } else {
      ll_info.convcode = 1;
    }
  } else {
    ll_info.convcode = -1 * convcode_ + 1;
  }

  ll_info.iters = as<arma::vec>(opt["iterations"])(0);
  
  if (ll_info.verbose) {
    Rcout << ll_info.LL << ' ';
    arma::vec& par(ll_info.min_par);
    for (uint_t i = 0; i < par.n_elem; i++) Rcout << par(i) << ' ';
    Rcout << std::endl;
  }
  
  return;
}


/*
 Fit `cor_phylo` model using R's `stats::optim`.
 
 Make sure this doesn't get run in parallel!

 */
void fit_cor_phylo_R(LogLikInfo& ll_info,
                     const double& rel_tol,
                     const int& max_iter,
                     const std::string& method,
                     const std::vector<double>& sann) {
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
  Rcpp::Function cor_phylo_LL_fxn = phyr_pkg["cor_phylo_LL"];
  
  Rcpp::List opt;

  NumericVector par0(ll_info.par0.begin(), ll_info.par0.end());
  
  if (method == "sann") {
    List sann_control;
    sann_control["maxit"] = sann[0];
    sann_control["temp"] = sann[1];
    sann_control["tmax"] = sann[2];
    sann_control["reltol"] = rel_tol;
    opt = optim(_["par"] = par0,
                _["fn"] = cor_phylo_LL_fxn,
                _["method"] = "SANN",
                _["control"] = sann_control,
                _["XX"] = ll_info.XX,
                _["UU"] = ll_info.UU,
                _["MM"] = ll_info.MM,
                _["Vphy"] = ll_info.Vphy,
                _["tau"] = ll_info.tau,
                _["REML"] = ll_info.REML,
                _["constrain_d"] = ll_info.constrain_d,
                _["lower_d"] = ll_info.lower_d,
                _["verbose"] = ll_info.verbose,
                _["rcond_threshold"] = ll_info.rcond_threshold);
    par0 = as<NumericVector>(opt["par"]);
  }
  
  List optim_control = List::create(_["maxit"] = max_iter, 
                                    _["reltol"] = rel_tol);
  
  opt = optim(_["par"] = par0,
              _["fn"] = cor_phylo_LL_fxn,
              _["method"] = "Nelder-Mead",
              _["control"] = optim_control,
              _["XX"] = ll_info.XX,
              _["UU"] = ll_info.UU,
              _["MM"] = ll_info.MM,
              _["Vphy"] = ll_info.Vphy,
              _["tau"] = ll_info.tau,
              _["REML"] = ll_info.REML,
              _["constrain_d"] = ll_info.constrain_d,
              _["lower_d"] = ll_info.lower_d,
              _["verbose"] = ll_info.verbose,
              _["rcond_threshold"] = ll_info.rcond_threshold);
  
  
  ll_info.min_par = as<arma::vec>(opt["par"]);
  
  ll_info.LL = as<double>(opt["value"]);
  ll_info.convcode = as<int>(opt["convergence"]);
  ll_info.iters = as<arma::vec>(opt["counts"])(0);
  
  if (ll_info.verbose) {
    Rcout << ll_info.LL << ' ';
    const arma::vec& par(ll_info.min_par);
    for (uint_t i = 0; i < par.n_elem; i++) Rcout << par(i) << ' ';
    Rcout << std::endl;
  }
  
  return;
  
}



/*
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 
 Other functions
 
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */

/*-
 Standardize matrices in place.
 
 Makes each column of the `X` matrix have mean of zero and standard deviation of 1.
 If `U` isn't empty, this function makes each column in each matrix have
 mean of zero and standard deviation of 1, unless all values are the same, in which
 case it keeps the standard deviation at zero.
 Divides each column of `M` by the original standard deviation of that column in 
 `X`.
 */
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
  
  for (uint_t i = 0; i < U.size(); i++) {
    arma::mat& Ui(U[i]);
    for (uint_t j = 0; j < Ui.n_cols; j++) {
      double sd = arma::stddev(Ui.col(j));
      Ui.col(j) -= arma::mean(Ui.col(j));
      if (sd > 0) Ui.col(j) /= sd;
    }
  }
  
  return;
}





/*
 Make an `LogLikInfo` object based on input matrices.
 The output `LogLikInfo` is used for model fitting.
 */
LogLikInfo::LogLikInfo(const arma::mat& X,
                 const std::vector<arma::mat>& U,
                 const arma::mat& M,
                 const arma::mat& Vphy_,
                 const bool& REML_,
                 const bool& no_corr_,
                 const bool& constrain_d_,
                 const double& lower_d_,
                 const bool& verbose_,
                 const double& rcond_threshold_) 
  : REML(REML_), no_corr(no_corr_), constrain_d(constrain_d_), lower_d(lower_d_),
    verbose(verbose_), rcond_threshold(rcond_threshold_), iters(0) {
  
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
  safe_chol(L, "model fitting");
  L = L.t();
  
  par0 = make_par(p, L, no_corr);
  min_par = par0;

}


/*
 Make an `LogLikInfo` object based on input matrices and another LogLikInfo object.
 
 The output `LogLikInfo` is used for model fitting.
 
 *Note:* This version is used for bootstrapping.
 It's different from the one above in that it doesn't re-normalize Vphy, UU, or tau.
 If you normalize Vphy and tau twice (which would happen if I used the previous
 version of this constructor), it can result in weird behavior.
 Notably, the bootstrap replicate will sometimes not converge, but when I output the
 same data and re-run cor_phylo on it, it'll converge.
 This is confusing, so I'm trying to avoid that.
 */
LogLikInfo::LogLikInfo(const arma::mat& X,
                 const std::vector<arma::mat>& U,
                 const arma::mat& M,
                 const LogLikInfo& other) 
  : UU(other.UU), Vphy(other.Vphy), tau(other.tau), REML(other.REML),
    no_corr(other.no_corr), constrain_d(other.constrain_d), lower_d(other.lower_d),
    verbose(other.verbose), rcond_threshold(other.rcond_threshold), iters(0) {

  uint_t p = X.n_cols;
  
  arma::mat Xs = X;
  std::vector<arma::mat> Us = U;
  arma::mat Ms = M;
  standardize_matrices(Xs, Us, Ms);
  
  
  XX = arma::reshape(Xs, Xs.n_elem, 1);
  MM = flex_pow(Ms, 2);
  MM.reshape(MM.n_elem, 1);
  
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
      }
    }
  }
  L = arma::cov(eps);
  safe_chol(L, "a bootstrap replicate");
  L = L.t();
  
  par0 = make_par(p, L, no_corr);
  min_par = par0;

}




inline void main_output(arma::mat& corrs, arma::mat& B, arma::mat& B_cov, arma::vec& d,
                        const LogLikInfo& ll_info,
                        const arma::mat& X, const std::vector<arma::mat>& U) {
  
  uint_t n = X.n_rows;
  uint_t p = X.n_cols;
  
  arma::mat L = make_L(ll_info.min_par, p);
  
  arma::mat R = L.t() * L;
  
  corrs = make_corrs(R);
  
  d = make_d(ll_info.min_par, p, ll_info.constrain_d, ll_info.lower_d);
  
  // OU transform
  arma::mat C = make_C(n, p, ll_info.tau, d, ll_info.Vphy, R);
  
  arma::mat V = make_V(C, ll_info.MM);
  
  arma::mat iV = arma::inv(V);
  
  arma::mat denom = ll_info.UU.t() * iV * ll_info.UU;
  
  arma::mat num = ll_info.UU.t() * iV * ll_info.XX;
  
  arma::vec B0 = arma::solve(denom, num);
  
  make_B_B_cov(B, B_cov, B0, iV, ll_info.UU, X, U);
  
  return;
}


/*
 Retrieve objects for output `cor_phylo` object.
 
 Returns a list containing output information, to later be coerced to a `cor_phylo`
   object by the `cor_phylo` function.
 */
List cp_get_output(const arma::mat& X,
                   const std::vector<arma::mat>& U,
                   const arma::mat& M,
                   const LogLikInfo& ll_info) {

  
  uint_t n = X.n_rows;
  uint_t p = X.n_cols;
  
  /*
   Get the main output from cor_phylo: correlations, coefficient estimates, 
   coefficient covariances, and phylogenetic signals.
   The objects constructed below correspond to this information.
   */
  arma::mat corrs;
  arma::mat B;
  arma::mat B_cov;
  arma::vec d;
  main_output(corrs, B, B_cov, d, ll_info, X, U);
  
  double logLik = -0.5 * std::log(2 * arma::datum::pi);
  if (ll_info.REML) {
    logLik *= (n * p - ll_info.UU.n_cols);
    arma::mat to_det = ll_info.XX.t() * ll_info.XX;
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
  
  std::vector<double> rcond_vals = return_rcond_vals(ll_info);
  
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
    _["rcond_vals"] = rcond_vals,
    _["min_par"] = ll_info.min_par,
    _["bootstrap"] = List::create()
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
//' @inheritParams cor_phylo
//' @param method the `method` input to `cor_phylo`.
//' 
//' @return a list containing output information, to later be coerced to a 
//'   `cor_phylo` object by the `cor_phylo` function.
//' @noRd
//' @name cor_phylo_cpp
//' 
//[[Rcpp::export]]
List cor_phylo_cpp(const arma::mat& X,
                   const std::vector<arma::mat>& U,
                   const arma::mat& M,
                   const arma::mat& Vphy,
                   const bool& REML,
                   const bool& constrain_d,
                   const double& lower_d,
                   const bool& verbose,
                   const double& rcond_threshold,
                   const double& rel_tol,
                   const int& max_iter,
                   const std::string& method,
                   const bool& no_corr,
                   const std::vector<double>& sann) {
  

  // LogLikInfo is C++ class to use for organizing info for optimizing
  LogLikInfo ll_info(X, U, M, Vphy, REML, no_corr, constrain_d, 
                     lower_d, verbose, rcond_threshold);

  /*
   Do the fitting.
   Methods "nelder-mead-r" and "sann" use R's `stats::optim`.
   Otherwise, use nlopt.
   */
  if (method == "nelder-mead-r" || method == "sann") {
    fit_cor_phylo_R(ll_info, rel_tol, max_iter, method, sann);
  } else {
    fit_cor_phylo_nlopt(ll_info, rel_tol, max_iter, method);
  }
  
  // Retrieve output from `ll_info` object and convert to list
  // Also do bootstrapping if desired
  List output = cp_get_output(X, U, M, ll_info);
  
  return output;
  
}

//' Inner function to do one bootstrapping rep.
//' 
//' The arguments before `X` are needed in this function but not in
//' `cor_phylo_cpp`, while the arguments starting at `X` are all those required
//' in `cor_phylo_cpp`.
//'
//' @noRd
//' 
//[[Rcpp::export]]
List one_boot_cpp(const arma::vec& rnd_vec,
                  const arma::vec& min_par,
                  const arma::mat& B,
                  const arma::vec& d,
                  const std::string& keep_boots,
                  const arma::mat& X,
                  const std::vector<arma::mat>& U,
                  const arma::mat& M,
                  const arma::mat& Vphy,
                  const bool& REML,
                  const bool& constrain_d,
                  const double& lower_d,
                  const bool& verbose,
                  const double& rcond_threshold,
                  const double& rel_tol,
                  const int& max_iter,
                  const std::string& method,
                  const bool& no_corr,
                  const std::vector<double>& sann) {
  
  uint_t n = X.n_rows;
  uint_t p = X.n_cols;
  
  LogLikInfo ll_info(X, U, M, Vphy, REML, no_corr, constrain_d, 
                     lower_d, verbose, rcond_threshold);
  ll_info.min_par = min_par;
  
  // Matrix to simulate error with phylogenetic signal
  arma::mat iD = make_iD(ll_info, p, d);
  
  // For predicted X values (i.e., without error)
  arma::mat pred_X = ll_info.UU.t();
  arma::rowvec tmp = B.col(0).t();
  pred_X = tmp * pred_X;
  pred_X = pred_X.t();
  pred_X.reshape(n, p);
  
  // Add error with phylogenetic signal (rnd_vec should be ~ N(0,1)):
  arma::mat rnd_X = iD * rnd_vec;
  rnd_X.reshape(n, p);
  for (uint_t i = 0; i < p; i++) {
    double sd_ = arma::stddev(X.col(i));
    rnd_X.col(i) *= sd_;
  }
  arma::mat new_X = pred_X + rnd_X;
  
  // Object to store info for new fit
  LogLikInfo ll_info_new(new_X, U, M, ll_info);
  
  /*
   Do the fitting.
   Methods "nelder-mead-r" and "sann" use R's `stats::optim`.
   Otherwise, use nlopt.
   */
  if (method == "nelder-mead-r" || method == "sann") {
    // Do the fitting:
    fit_cor_phylo_R(ll_info_new, rel_tol, max_iter, method, sann);
  } else {
    // Same thing for nlopt version
    fit_cor_phylo_nlopt(ll_info_new, rel_tol, max_iter, method);
  }
  
  // Primary output objects:
  arma::mat new_corrs;
  arma::mat new_B;
  arma::mat new_B_cov;
  arma::vec new_d;
  main_output(new_corrs, new_B, new_B_cov, new_d, ll_info_new, new_X, U);
  
  List boot_list = List::create(_["corrs"] = new_corrs, 
                                _["d"] = new_d,
                                _["B0"] = new_B.col(0), 
                                _["B_cov"] = new_B_cov,
                                _["convcodes"] = ll_info_new.convcode);
  
  // Determine whether convergence failed, and save simulated data if specified
  bool failed = ll_info_new.convcode != 0;
  if (keep_boots == "all" || (keep_boots == "fail" && failed)) {
    boot_list["mats"] = new_X;
  } else if (keep_boots == "fail" && !failed) {
    // Add filler matrix if only keeping failed reps and this one succeeded
    boot_list["mats"] = arma::mat(0,0);
  }
  
  return boot_list;
  
}


// Organize bootstrap output.
//[[Rcpp::export]]
List organize_boots(List orig_list) {
  
  uint_t n_boots = orig_list.size();
  
  // Get some basic info from the first bootstrap
  List front_boot(orig_list[0]);
  bool has_mats = front_boot.containsElementNamed("mats");
  uint_t p = as<NumericMatrix>(front_boot["corrs"]).nrow();
  uint_t B_rows = as<NumericVector>(front_boot["B0"]).size();
  
  arma::cube corrs(p, p, n_boots);
  arma::mat d(p, n_boots);
  arma::mat B0(B_rows, n_boots);
  arma::cube B_cov(B_rows, B_rows, n_boots);
  arma::Col<int> convcodes(n_boots);
  List mats;
  if (has_mats) mats = List(n_boots);
  
  arma::mat corrs_i;
  arma::vec B0_i;
  arma::mat B_cov_i;
  arma::vec d_i;
  NumericMatrix mats_i;
  
  for (uint_t i = 0; i < n_boots; i++) {
    
    List boot_i(orig_list[i]);
    
    corrs_i = as<arma::mat>(boot_i["corrs"]);
    B0_i = as<arma::vec>(boot_i["B0"]);
    B_cov_i = as<arma::mat>(boot_i["B_cov"]);
    d_i = as<arma::vec>(boot_i["d"]);
    
    corrs.slice(i) = corrs_i;
    B0.col(i) = B0_i;
    B_cov.slice(i) = B_cov_i;
    d.col(i) = d_i;
    convcodes(i) = as<int>(boot_i["convcodes"]);
    
    if (has_mats) {
      mats_i = as<NumericMatrix>(boot_i["mats"]);
      mats(i) = mats_i;
    }
    
  }
  
  List new_list = List::create(_["corrs"] = corrs,
                               _["d"] = d,
                               _["B0"] = B0, 
                               _["B_cov"] = B_cov,
                               _["convcodes"] = convcodes,
                               _["mats"] = mats);
  
  return new_list;
}

