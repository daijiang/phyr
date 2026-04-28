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
arma::vec pglmm_gaussian_predict(const arma::mat& iV,
                                 const arma::mat& H){
  // Conditional mean E[H_i | H_{-i}] via precision matrix identity:
  //   V[i,-i] * V[-i,-i]^{-1} = -iV[i,-i] / iV[ii]
  // so h(i) = H[i] - (iV * H)[i] / iV[ii]
  // Cost: O(n^2) matvec vs the original O(n^4) loop of n matrix inversions.
  arma::vec Hv      = arma::vectorise(H);   // handles n×1 mat or vec input
  arma::vec iV_H    = iV * Hv;              // O(n^2) matvec
  arma::vec diag_iV = iV.diag();            // O(n)
  return Hv - iV_H / diag_iV;              // O(n) element-wise
}

// [[Rcpp::export]]
double pglmm_gaussian_LL_cpp(NumericVector par,
                           const arma::mat& X, const arma::vec& Y,
                           const arma::sp_mat& Zt, const arma::sp_mat& St,
                           const List& nested,
                           bool REML, bool verbose){
  int n = X.n_rows;
  int p = X.n_cols;
  int q_nonNested = St.n_rows;
  arma::sp_mat Ut;
  arma::sp_mat U;
  if(q_nonNested > 0){
    IntegerVector idx = seq_len(q_nonNested) - 1; // c++ starts with 0
    NumericVector sr0 = par[idx];
    rowvec sr = real(as<rowvec>(sr0));
    arma::mat iC0 = sr * St;
    arma::vec iC1 = vectorise(iC0, 0); // extract by columns
    arma::sp_mat iC = sp_mat(diagmat(iC1));
    Ut = iC * Zt;
    U = trans(Ut);
  }
  int q_Nested = nested.size();

  NumericVector sn; // pre-declare out of if{}
  if (q_Nested > 0) {
    IntegerVector idx2 = wrap(seq(q_nonNested, q_nonNested + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    rowvec sn1 = real(as<rowvec>(sn0));
    sn = as<NumericVector>(wrap(sn1));
  }

  // Variables filled in each branch below
  double logdetV = 0.0;
  double HiVH   = 0.0;
  arma::mat denom;  // X' iV X (p x p), needed for REML

  if (q_Nested == 0) {
    // ---- non-nested only: A = I, matrix-free Woodbury (OPT-G0) ----
    // iV = I - U*(I+Ut*U)^{-1}*Ut; never form n×n iV.
    // iV*x = x - U*(M*(Ut*x)), O(n*Q) vs old O(n^2).
    arma::mat Ut_d = arma::mat(Ut);   // Q × n
    arma::mat U_d  = Ut_d.t();        // n × Q
    int Q = (int)Ut_d.n_rows;

    arma::mat IpUtU = arma::eye(Q, Q) + Ut_d * U_d;   // Q × Q
    arma::mat M     = arma::inv_sympd(IpUtU);

    // logdetV = log|I + Ut*U|  (logdet(I_n)=0)
    double signV; arma::log_det(logdetV, signV, IpUtU);
    if (!std::isfinite(logdetV)) {
      arma::mat lgm = arma::chol(IpUtU);
      logdetV = 2.0 * arma::accu(arma::log(lgm.diag()));
    }

    // iV*x helpers via Woodbury (no n×n matrix)
    auto iV_vec = [&](const arma::vec& x) -> arma::vec {
      return x - U_d * (M * (Ut_d * x));
    };
    auto iV_mat = [&](const arma::mat& A2) -> arma::mat {
      return A2 - U_d * (M * (Ut_d * A2));
    };

    arma::mat iVX  = iV_mat(X);
    denom          = X.t() * iVX;
    arma::mat num  = X.t() * iV_vec(Y);
    arma::mat B    = arma::solve(denom, num);
    arma::vec H    = Y - X * B;
    arma::vec iVH  = iV_vec(H);
    HiVH           = arma::dot(H, iVH);

  } else {
    // ---- nested terms present ----
    // Build A = I + sum_j sn_j^2 * nested[j]
    arma::sp_mat A = sp_mat(n, n); A.eye();
    if (q_Nested == 1){
      double snj = pow(sn[0], 2);
      sp_mat nj = nested[0];
      A = A + snj * nj;
    } else {
      for (int j = 0; j < q_Nested; j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
    arma::mat A1(A);

    // OPT1+OPT2: Cholesky once for logdetA; then use triangular solves to
    // compute all iA*x products without ever materialising the n×n matrix iA.
    // Savings vs old LU+log_det(iV): eliminates ~4/3 n^3 FLOPS per LL call.
    arma::mat chol_A;
    bool chol_ok = arma::chol(chol_A, A1);  // upper R s.t. R'R = A

    if (chol_ok) {
      double logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
      // Pre-transpose R once so trimatl() reuses it without repeated copies.
      arma::mat Rlo = chol_A.t();  // lower-triangular L s.t. L*L' = A (R')

      // Solve A*Z = RHS using stored factors:
      //   forward:  Rlo * W = RHS  (lower-tri solve)
      //   backward: chol_A * Z = W (upper-tri solve)
      arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
      iA_X           = arma::solve(arma::trimatu(chol_A), iA_X);

      arma::mat iA_Y = arma::solve(arma::trimatl(Rlo), arma::mat(Y));
      iA_Y           = arma::solve(arma::trimatu(chol_A), iA_Y);

      if (q_nonNested > 0) {
        // Woodbury: iV = iA - iA U M Ut iA,  M = (I + Ut iA U)^{-1}
        arma::mat U_dense(U);
        arma::mat iA_U = arma::solve(arma::trimatl(Rlo), U_dense);
        iA_U           = arma::solve(arma::trimatu(chol_A), iA_U);

        arma::mat Ut_dense(Ut);
        arma::mat Ut_iA_U = Ut_dense * iA_U;                   // Q x Q (Q = sum of levels)
        int Q = Ut.n_rows;
        arma::mat Ishort_Ut_iA_U = arma::eye(Q, Q) + Ut_iA_U;
        arma::mat M = arma::inv_sympd(Ishort_Ut_iA_U);         // Q x Q

        // X' iV X = iA_X' X - (Ut iA_X)' M (Ut iA_X)
        arma::mat Ut_iA_X = Ut_dense * iA_X;                   // q x p
        arma::vec Ut_iA_Y = Ut_dense * iA_Y.col(0);            // q x 1

        denom          = iA_X.t() * X  - Ut_iA_X.t() * M * Ut_iA_X;
        arma::vec num  = iA_X.t() * Y  - Ut_iA_X.t() * M * Ut_iA_Y;
        arma::vec B    = arma::solve(denom, num);
        arma::vec H    = Y - X * B;

        // H' iV H = dot(H, iA_H) - (Ut iA_H)' M (Ut iA_H)
        arma::mat iA_H = arma::solve(arma::trimatl(Rlo), arma::mat(H));
        iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
        arma::vec Ut_iA_H = Ut_dense * iA_H.col(0);
        HiVH = arma::dot(H, iA_H.col(0)) -
               arma::as_scalar(Ut_iA_H.t() * M * Ut_iA_H);

        // logdetV via matrix-det lemma: log|V| = log|A| + log|I + Ut A^{-1} U|
        double signV;
        arma::log_det(logdetV, signV, Ishort_Ut_iA_U);
        logdetV += logdetA;

      } else {
        // q_nonNested == 0: iV = iA
        denom          = iA_X.t() * X;
        arma::vec num  = iA_X.t() * Y;
        arma::vec B    = arma::solve(denom, num);
        arma::vec H    = Y - X * B;

        arma::mat iA_H = arma::solve(arma::trimatl(Rlo), arma::mat(H));
        iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
        HiVH    = arma::dot(H, iA_H.col(0));
        logdetV = logdetA;
      }

    } else {
      // Cholesky failed (A not PD for these params): fallback to LU
      arma::mat iA1 = arma::inv(A1);
      double logdetA_fb, sgn;
      arma::log_det(logdetA_fb, sgn, A1);
      arma::sp_mat iA = sp_mat(iA1);

      arma::sp_mat iV0;
      arma::mat Ishort_Ut_iA_U_fb;
      if (q_nonNested > 0) {
        arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
        arma::sp_mat Ut_iA_U = Ut * iA * U;
        Ishort_Ut_iA_U_fb = mat(Ishort + Ut_iA_U);
        arma::mat i_Ishort = arma::inv_sympd(Ishort_Ut_iA_U_fb);
        iV0 = iA - iA * U * sp_mat(i_Ishort) * Ut * iA;
      } else {
        iV0 = iA;
      }
      arma::mat iV(iV0);
      denom         = trans(X) * iV * X;
      arma::mat num = trans(X) * iV * Y;
      arma::mat B   = solve(denom, num);
      arma::vec H   = Y - X * B;
      HiVH          = as_scalar(trans(H) * iV * H);

      if (q_nonNested > 0) {
        double signV;
        log_det(logdetV, signV, Ishort_Ut_iA_U_fb);
        logdetV += logdetA_fb;
      } else {
        logdetV = logdetA_fb;
      }
    }
  }

  // ---- Compute LL from branch-agnostic quantities ----
  double LL;
  if(REML){
    double s2_conc = HiVH / (n - p);
    double logdetL, signL;
    arma::log_det(logdetL, signL, denom);
    LL = 0.5 * ((n - p) * std::log(s2_conc) + logdetV + (n - p) + logdetL);
  } else {
    double s2_conc = HiVH / n;
    LL = 0.5 * (n * std::log(s2_conc) + logdetV + n);
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
  arma::sp_mat Ut;
  arma::sp_mat U;
  arma::rowvec sr;
  if(q_nonNested > 0){
    IntegerVector idx = seq_len(q_nonNested) - 1; // c++ starts with 0
    // uvec idx_uvec = as<uvec>(idx);
    NumericVector sr0 = par[idx];
    sr = real(as<rowvec>(sr0));
    arma::mat iC0 = sr * St;
    arma::vec iC1 = vectorise(iC0, 0); // extract by columns
    arma::sp_mat iC = sp_mat(diagmat(iC1));
    Ut = iC * Zt;
    U = trans(Ut);
  }
  
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
  double logdetA = 0.0;  // in scope for log-det block below
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
      for (int j = 0; j < q_Nested; j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
    arma::mat A1(A);
    // OPT1: single Cholesky for both logdetA and iA
    arma::mat chol_A;
    bool chol_ok = arma::chol(chol_A, A1);
    arma::sp_mat iA;
    if (chol_ok) {
      logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
      arma::mat Ri = arma::inv(trimatu(chol_A));
      iA = sp_mat(Ri * Ri.t());
    } else {
      arma::mat iA1 = arma::inv(A1);
      double sgn;
      arma::log_det(logdetA, sgn, A1);
      iA = sp_mat(iA1);
    }
    if(q_nonNested > 0){
      arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
      arma::sp_mat Ut_iA_U = Ut * iA * U;
      Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
      arma::mat i_Ishort_Ut_iA_U = arma::inv_sympd(Ishort_Ut_iA_U);
      iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
    } else {
      iV0 = iA;
    }
  }

  arma::mat iV(iV0); // convert to dense matrix
  arma::mat denom = trans(X) * iV * X;
  arma::mat num = trans(X) * iV * Y;
  arma::mat B = solve(denom, num);
  arma::vec H = Y - X * B;

  // logdetV is computed but not returned; kept for consistency
  double logdetV;
  double signV;
  if (q_Nested == 0) {
    log_det(logdetV, signV, Ishort_Ut_iA_U);
    NumericVector logdetV1 = NumericVector::create(logdetV);
    if(any(is_infinite(logdetV1))){
      arma::mat lgm = chol(Ishort_Ut_iA_U);
      logdetV = 2 * sum(log(lgm.diag()));
    }
  } else {
    // OPT1 (cont.): matrix-determinant lemma
    if (q_nonNested > 0) {
      log_det(logdetV, signV, Ishort_Ut_iA_U);
      logdetV += logdetA;
    } else {
      logdetV = logdetA;
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
  Rcpp::checkUserInterrupt();
  // start optimization
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"]; 
  Rcpp::Environment nloptr_pkg = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptr_pkg["nloptr"];
  Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
  Rcpp::Function pglmm_gaussian_LL_cpp_fxn = phyr_pkg["pglmm_gaussian_LL_cpp"];
  
  Rcpp::List opt;
  if(optimizer == "Nelder-Mead"){
    if(q > 1){
      opt = optim(_["par"]    = par,
                  _["fn"]     = pglmm_gaussian_LL_cpp_fxn,
                  _["X"] = X, _["Y"] = Y, _["Zt"] = Zt,
                  _["St"] = St, _["nested"] = nested,
                  _["REML"] = REML, _["verbose"] = verbose,
                  _["method"] = "Nelder-Mead",
                  _["control"] = List::create(_["maxit"] = maxit, _["reltol"] = reltol));
    } else {
      opt = optim(_["par"]    = par,
                  _["fn"]     = pglmm_gaussian_LL_cpp_fxn,
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
                     _["eval_f"] = pglmm_gaussian_LL_cpp_fxn,
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

/*** R
# pglmm_gaussian_predict(x$iV, x$H)
# pglmm_gaussian_internal_cpp(par = s, X, Y, Zt = as(matrix(0, 0, 0), "dgTMatrix"), 
#                             St = as(matrix(0, 0, 0), "dgTMatrix"), nested, REML, 
#                             verbose, optimizer, maxit, 
#                             reltol, q, n, p, pi)
# res = pglmm_gaussian_internal_cpp(par = s, X, Y, Zt, St, nested, REML, 
#                             verbose, optimizer, maxit, 
#                             reltol, q, n, p, pi)
*/