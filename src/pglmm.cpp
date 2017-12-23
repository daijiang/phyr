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
    Rcout << LL << " " << par << " ";
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
List plmm_binary_iV_logdetV_cpp(NumericVector par, arma::vec mu,
                                const arma::sp_mat& Zt, const arma::sp_mat& St, 
                                const List& nested, bool logdet){
  
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
  double logdetV;
  double signV;
  double logdetiA;
  double signiA;
  if (q_Nested == 0){
    arma::vec pq = mu % (1 - mu);
    arma::sp_mat iA = sp_mat(diagmat(pq));
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * iA * U;
    // Woodbury identity
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
    arma::mat iV(iV0); // convert to dense matrix
    if(logdet){
      log_det(logdetV, signV, Ishort_Ut_iA_U); 
      log_det(logdetiA, signiA, mat(iA)); 
      logdetV = logdetV - logdetiA;
      NumericVector logdetV1 = NumericVector::create(logdetV);
      if(any(is_infinite(logdetV1))){
        arma::mat lgm = chol(Ishort_Ut_iA_U);
        logdetV = 2 * sum(log(lgm.diag())) - logdetiA;
      }
    }
  } else {
    arma::vec pq = 1 / (mu % (1 - mu));
    arma::sp_mat A = sp_mat(diagmat(pq));
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
    arma::mat iV(iV0); // convert to dense matrix
    if(logdet){
      log_det(logdetV, signV, iV); 
      logdetV = -1 * logdetV;
      NumericVector logdetV1 = NumericVector::create(logdetV);
      if(any(is_infinite(logdetV1))){
        arma::mat lgm = chol(iV);
        logdetV = -2 * sum(log(lgm.diag()));
      }
    }
  }

  
  if(logdet){
    return List::create(
      _["iV"] = iV0,
      _["logdetV"] = logdetV
    );
  } else {
    return List::create(_["iV"] = iV0);
  }
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# par = c(1.002103, 1.002103, 1.002103, 1.002103)
# 
# x = pglmm_gaussian_LL_calc_cpp(par = par, X = X, Y = Y, Zt = Zt, St = St, nested = nested, REML = T)
# x
# expect_equal(x, as.numeric(LL))
# expect_equal(x$Ut, Ut)
# expect_equal(x$U, U)
# expect_equivalent(x$sn[1,1], sn)
# expect_equivalent(x$iV, as.matrix(iV))
# expect_equivalent(x$B, B)
# expect_equivalent(x$H, H)
# expect_equivalent(x$logdetV, logdetV)
x = plmm_binary_iV_logdetV_cpp(par = ss, Zt = Zt, St = St, mu = mu, nested = nested, logdet = T)
expect_equal(as.matrix(x$iV), as.matrix(y$iV))
expect_equal(x$logdetV, y$logdetV)
*/
