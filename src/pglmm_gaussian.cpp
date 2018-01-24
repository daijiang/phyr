// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double pglmm_gaussian_LL_cpp(NumericVector par, 
                           const arma::mat& X, const arma::vec& Y, 
                           const arma::sp_mat& Zt, const arma::sp_mat& St, 
                           const List& nested, 
                           bool REML, bool verbose){
  int n = X.n_rows;
  int p = X.n_cols;
  int q_nonNested = St.n_rows;
  IntegerVector idx = seq_len(q_nonNested) - 1; // c++ starts with 0
  // uvec idx_uvec = as<uvec>(idx);
  NumericVector sr0 = par[idx];
  rowvec sr = real(as<rowvec>(sr0));
  arma::mat iC0 = sr * St;
  arma::vec iC1 = vectorise(iC0, 0); // extract by columns
  arma::sp_mat iC = sp_mat(diagmat(iC1));
  arma::sp_mat Ut = iC * Zt;
  arma::sp_mat U = trans(Ut);
  int q_Nested = nested.size();

  NumericVector sn; // pre-declare out of if{}
  if (q_Nested > 0) {
    IntegerVector idx2 = wrap(seq(q_nonNested, q_nonNested + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    rowvec sn1 = real(as<rowvec>(sn0));
    sn = as<NumericVector>(wrap(sn1)); // no need to declare type again
  } 
  
  arma::sp_mat iV0;
  arma::mat Ishort_Ut_iA_U;
  if (q_Nested == 0){
    arma::sp_mat iA = sp_mat(n, n); iA.eye();
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * U;
    // Woodbury identity
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - U * sp_mat(i_Ishort_Ut_iA_U) * Ut;
  } else {
    arma::sp_mat A = sp_mat(n, n); A.eye();
    if (q_Nested == 1){
      double snj = pow(sn[0], 2);
      sp_mat nj = nested[0];
      A = A + snj * nj;
    } else {
      for (int j = 0; j < (q_Nested - 1); j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
    arma::mat A1(A);
    arma::sp_mat iA = sp_mat(inv(A1));
    // Rcout << iA << " " ;
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * iA * U;
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
  }
  
  arma::mat iV(iV0); // convert to dense matrix
  arma::mat denom = trans(X) * iV * X;
  arma::mat num = trans(X) * iV * Y;
  arma::mat B = solve(denom, num);
  arma::vec H = Y - X * B;
  
  double logdetV;
  double signV;
  if (q_Nested == 0) {
    // Sylvester identity
    log_det(logdetV, signV, Ishort_Ut_iA_U); 
    NumericVector logdetV1 = NumericVector::create(logdetV);
    if(any(is_infinite(logdetV1))){
      arma::mat lgm = chol(Ishort_Ut_iA_U);
      logdetV = 2 * sum(log(lgm.diag()));
    }
  } else {
    log_det(logdetV, signV, iV);
    logdetV = -1 * logdetV;
    NumericVector logdetV1 = NumericVector::create(logdetV);
    if(any(is_infinite(logdetV1))){
      arma::mat lgm = chol(iV);
      logdetV = -2 * sum(log(lgm.diag()));
    }
  }
  
  double LL;
  if(REML){
    double s2_conc = as_scalar(trans(H) * iV * H) / (n - p);
    double logdetL;
    double signL;
    log_det(logdetL, signL, denom);
    LL = 0.5 * ((n - p) * log(s2_conc) + logdetV + (n - p) + logdetL);
  } else {
    double s2_conc = as_scalar(trans(H) * iV * H) / n;
    LL = 0.5 * (n * log(s2_conc) + logdetV + n);
  }
  
  if(verbose){
    Rcout << LL << " " << par << std::endl;
  }
    
  return LL;
}

// [[Rcpp::export]]
List pglmm_gaussian_LL_calc_cpp(NumericVector par, 
                                const arma::mat& X, const arma::vec& Y, 
                                const arma::sp_mat& Zt, const arma::sp_mat& St, 
                                const List& nested, bool REML){
  int n = X.n_rows;
  int p = X.n_cols;
  int q_nonNested = St.n_rows;
  IntegerVector idx = seq_len(q_nonNested) - 1; // c++ starts with 0
  // uvec idx_uvec = as<uvec>(idx);
  NumericVector sr0 = par[idx];
  rowvec sr = real(as<rowvec>(sr0));
  arma::mat iC0 = sr * St;
  arma::vec iC1 = vectorise(iC0, 0); // extract by columns
  arma::sp_mat iC = sp_mat(diagmat(iC1));
  arma::sp_mat Ut = iC * Zt;
  arma::sp_mat U = trans(Ut);
  int q_Nested = nested.size();
  
  NumericVector sn; // pre-declare out of if{}
  if (q_Nested > 0) {
    IntegerVector idx2 = wrap(seq(q_nonNested, q_nonNested + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    rowvec sn1 = real(as<rowvec>(sn0));
    sn = as<NumericVector>(wrap(sn1)); // no need to declare type again
  } 
  
  arma::sp_mat iV0;
  arma::mat Ishort_Ut_iA_U;
  if (q_Nested == 0){
    arma::sp_mat iA = sp_mat(n, n); iA.eye();
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * U;
    // Woodbury identity
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - U * sp_mat(i_Ishort_Ut_iA_U) * Ut;
  } else {
    arma::sp_mat A = sp_mat(n, n); A.eye();
    if (q_Nested == 1){
      double snj = pow(sn[0], 2);
      sp_mat nj = nested[0];
      A = A + snj * nj;
    } else {
      for (int j = 0; j < (q_Nested - 1); j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
    arma::mat A1(A);
    arma::sp_mat iA = sp_mat(inv(A1));
    // Rcout << iA << " " ;
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * iA * U;
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
  }
  
  arma::mat iV(iV0); // convert to dense matrix
  arma::mat denom = trans(X) * iV * X;
  arma::mat num = trans(X) * iV * Y;
  arma::mat B = solve(denom, num);
  arma::vec H = Y - X * B;
  
  double logdetV;
  double signV;
  if (q_Nested == 0) {
    // Sylvester identity
    log_det(logdetV, signV, Ishort_Ut_iA_U); 
    NumericVector logdetV1 = NumericVector::create(logdetV);
    if(any(is_infinite(logdetV1))){
      arma::mat lgm = chol(Ishort_Ut_iA_U);
      logdetV = 2 * sum(log(lgm.diag()));
    }
  } else {
    log_det(logdetV, signV, iV);
    logdetV = -1 * logdetV;
    NumericVector logdetV1 = NumericVector::create(logdetV);
    if(any(is_infinite(logdetV1))){
      arma::mat lgm = chol(iV);
      logdetV = -2 * sum(log(lgm.diag()));
    }
  }
  
  double s2resid;
  if(REML){
    s2resid = as_scalar(trans(H) * iV * H) / (n - p);
  } else {
    s2resid = as_scalar(trans(H) * iV * H) / n;
  }
  
  iV = iV/s2resid;
  rowvec s2r = s2resid * pow(sr, 2);
  NumericVector s2n = s2resid * pow(sn, 2);
  arma::mat B_cov = inv(trans(X) * iV * X);
  arma::vec B_se = sqrt(B_cov.diag());
  
  return List::create(
    _["B"] = B,
    _["B.se"] = B_se,
    _["B.cov"] = B_cov,
    _["sr"] = sr,
    _["sn"] = sn,
    _["s2n"] = s2n,
    _["s2r"] = s2r,
    _["s2resid"] = s2resid,
    _["iV"] = iV,
    _["H"] = H
  );
}

// [[Rcpp::export]]
Rcpp::List pglmm_gaussian_internal_cpp(NumericVector par, 
                                       const arma::mat& X, const arma::vec& Y, 
                                       const arma::sp_mat& Zt, const arma::sp_mat& St, 
                                       const List& nested, bool REML, bool verbose,
                                       std::string optimizer, int maxit, double reltol,
                                       int q, int n, int p, const double Pi
                                       ){
  // start optimization
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"]; 
  Rcpp::Environment nloptr_pkg = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptr_pkg["nloptr"];
  
  Rcpp::List opt;
  if(optimizer == "Nelder-Mead"){
    if(q > 1){
      opt = optim(_["par"]    = par,
                  _["fn"]     = Rcpp::InternalFunction(&pglmm_gaussian_LL_cpp),
                  _["X"] = X, _["Y"] = Y, _["Zt"] = Zt,
                  _["St"] = St, _["nested"] = nested,
                  _["REML"] = REML, _["verbose"] = verbose,
                  _["method"] = "Nelder-Mead",
                  _["control"] = List::create(_["maxit"] = maxit, _["reltol"] = reltol));
    } else {
      opt = optim(_["par"]    = par,
                  _["fn"]     = Rcpp::InternalFunction(&pglmm_gaussian_LL_cpp),
                  _["X"] = X, _["Y"] = Y, _["Zt"] = Zt,
                  _["St"] = St, _["nested"] = nested,
                  _["REML"] = REML, _["verbose"] = verbose,
                  _["method"] = "L-BFGS-B",
                  _["control"] = List::create(_["maxit"] = maxit));
    }
  } else {
    std::string nlopt_algor;
    if (optimizer == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA";
    if (optimizer == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD";
    if (optimizer == "subplex") nlopt_algor = "NLOPT_LN_SBPLX";
    List opts = List::create(_["algorithm"] = nlopt_algor,
                             _["ftol_rel"] = reltol, _["ftol_abs"] = reltol,
                             _["xtol_rel"] = 0.0001,
                             _["maxeval"] = maxit);
    List S0 = nloptr(_["x0"] = par,
                     _["eval_f"] = Rcpp::InternalFunction(&pglmm_gaussian_LL_cpp),
                     _["opts"] = opts, _["X"] = X, _["Y"] = Y, _["Zt"] = Zt,
                     _["St"] = St, _["nested"] = nested,
                     _["REML"] = REML, _["verbose"] = verbose);
    opt = List::create(_["par"] = S0["solution"], _["value"] = S0["objective"],
                      _["counts"] = S0["iterations"], _["convergence"] = S0["status"],
                      _["message"] = S0["message"]);
  }
  // end of optimization
  
  arma::vec par_opt0 = abs(real(as<arma::vec>(opt["par"])));
  NumericVector par_opt = wrap(par_opt0);
  double LL = as_scalar(as<double>(opt["value"]));
  int convcode = as<int>(opt["convergence"]);
  arma::vec niter = as<arma::vec>(opt["counts"]);
  
  // calculate coef
  List out = pglmm_gaussian_LL_calc_cpp(par_opt, X, Y, Zt, St, nested, REML);
  double logLik, detx, signx;
  if(REML){
    log_det(detx, signx, trans(X) * X);
    logLik = -0.5 * (n - p) * log(2 * Pi) + 0.5 * detx - LL;
  } else {
    logLik = -0.5 * n * log(2 * Pi) - LL;
  }
  
  // return results
  return List::create(_["out"] = out, _["logLik"] = logLik,
                      _["convcode"] = convcode, _["niter"] = niter);
}

