// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "pglmm_utils.h"   // detect_n_sp, verify_block_structure, block_chol

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

  // OPT-C: detect block-diagonal structure once (geometry fixed per fit).
  // Block-Chol is only used for the nested-only path (q_nonNested == 0).
  int  n_sp_g   = 0, n_site_g = 0;
  bool use_blk_g = false;
  if (q_Nested > 0 && q_nonNested == 0) {
    arma::sp_mat nj0_sp = nested[0];
    n_sp_g   = detect_n_sp(nj0_sp);
    n_site_g = (n_sp_g > 0) ? n / n_sp_g : 0;
    use_blk_g = (n_sp_g > 0 && n_sp_g * n_site_g == n &&
                 verify_block_structure(nj0_sp, n_sp_g));
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

    if (q_nonNested == 0 && use_blk_g) {
      // OPT-C: block-diagonal Cholesky — O(n_site × n_sp³) vs O(n³).
      // Restricted to nested-only (q_nonNested == 0): the Woodbury step amplifies
      // tiny iA differences between block and full-dense paths when q_nonNested > 0.

      // Precompute per-block nested contribution (sn fixed for this LL call)
      std::vector<arma::mat> nj_fixed(n_site_g,
                                      arma::mat(n_sp_g, n_sp_g, arma::fill::zeros));
      for (int j = 0; j < q_Nested; ++j) {
        double snj2 = pow((double)sn[j], 2);
        arma::sp_mat njsp = nested[j];
        for (int k = 0; k < n_site_g; ++k) {
          int s = k * n_sp_g;
          nj_fixed[k] += snj2 *
            arma::mat(njsp(arma::span(s, s + n_sp_g - 1),
                           arma::span(s, s + n_sp_g - 1)));
        }
      }

      // pq = ones for Gaussian (A_k = I_{n_sp} + nj_fixed[k])
      arma::vec pq_ones(n, arma::fill::ones);
      std::vector<arma::mat> chols_blk(n_site_g);
      double logdetA = block_chol(chols_blk, pq_ones, nj_fixed, n_sp_g, n_site_g);

      // Block triangular solve helpers (per-block forward+back solve)
      auto blk_solve_v = [&](const arma::vec& rhs) -> arma::vec {
        arma::vec out(n);
        for (int k = 0; k < n_site_g; ++k) {
          int s = k * n_sp_g;
          arma::mat Rlo_k = chols_blk[k].t();
          arma::vec tmp = arma::solve(arma::trimatl(Rlo_k),
                                      rhs.subvec(s, s + n_sp_g - 1));
          out.subvec(s, s + n_sp_g - 1) = arma::solve(arma::trimatu(chols_blk[k]), tmp);
        }
        return out;
      };
      auto blk_solve_m = [&](const arma::mat& RHS) -> arma::mat {
        arma::mat out(RHS.n_rows, RHS.n_cols);
        for (int k = 0; k < n_site_g; ++k) {
          int s = k * n_sp_g;
          arma::mat Rlo_k = chols_blk[k].t();
          arma::mat tmp = arma::solve(arma::trimatl(Rlo_k),
                                      RHS.rows(s, s + n_sp_g - 1));
          out.rows(s, s + n_sp_g - 1) = arma::solve(arma::trimatu(chols_blk[k]), tmp);
        }
        return out;
      };

      arma::mat iA_X  = blk_solve_m(X);
      arma::vec iA_Y  = blk_solve_v(Y);
      denom           = iA_X.t() * X;
      arma::vec num   = iA_X.t() * Y;
      arma::vec B     = arma::solve(denom, num);
      arma::vec H     = Y - X * B;
      arma::vec iA_H  = blk_solve_v(H);
      HiVH    = arma::dot(H, iA_H);
      logdetV = logdetA;

    } else {
      // Full dense path: Build A = I + sum_j sn_j^2 * nested[j]
      arma::sp_mat A = sp_mat(n, n); A.eye();
      for (int j = 0; j < q_Nested; j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
      arma::mat A1(A);

      // OPT1+OPT2: Cholesky once for logdetA; then use triangular solves to
      // compute all iA*x products without ever materialising the n×n matrix iA.
      arma::mat chol_A;
      bool chol_ok = arma::chol(chol_A, A1);  // upper R s.t. R'R = A

      if (chol_ok) {
        double logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
        arma::mat Rlo = chol_A.t();  // lower L s.t. L*L' = A

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
          arma::mat Ut_iA_U = Ut_dense * iA_U;
          int Q = Ut.n_rows;
          arma::mat Ishort_Ut_iA_U = arma::eye(Q, Q) + Ut_iA_U;
          arma::mat M = arma::inv_sympd(Ishort_Ut_iA_U);

          arma::mat Ut_iA_X = Ut_dense * iA_X;
          arma::vec Ut_iA_Y = Ut_dense * iA_Y.col(0);

          denom          = iA_X.t() * X  - Ut_iA_X.t() * M * Ut_iA_X;
          arma::vec num  = iA_X.t() * Y  - Ut_iA_X.t() * M * Ut_iA_Y;
          arma::vec B    = arma::solve(denom, num);
          arma::vec H    = Y - X * B;

          arma::mat iA_H = arma::solve(arma::trimatl(Rlo), arma::mat(H));
          iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
          arma::vec Ut_iA_H = Ut_dense * iA_H.col(0);
          HiVH = arma::dot(H, iA_H.col(0)) -
                 arma::as_scalar(Ut_iA_H.t() * M * Ut_iA_H);

          // logdetV via matrix-det lemma
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
        // Cholesky failed: fallback to LU
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
    } // end full-dense branch
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

  // ---- build Ut / U ----
  arma::sp_mat Ut, U;
  arma::rowvec sr;
  if(q_nonNested > 0){
    IntegerVector idx = seq_len(q_nonNested) - 1;
    NumericVector sr0 = par[idx];
    sr = real(as<rowvec>(sr0));
    arma::vec iC1 = vectorise(sr * St, 0);
    Ut = sp_mat(diagmat(iC1)) * Zt;
    U  = trans(Ut);
  }

  int q_Nested = nested.size();
  NumericVector sn;
  if(q_Nested > 0){
    IntegerVector idx2 = wrap(seq(q_nonNested, q_nonNested + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    sn = as<NumericVector>(wrap(real(as<rowvec>(sn0))));
  }

  // Branch-agnostic outputs
  arma::vec B, H;
  double HiVH = 0.0, s2resid = 0.0;
  arma::mat denom_scaled;   // X'*(iV/s2resid)*X  for B_cov
  arma::mat iV_scaled;      // full dense iV/s2resid returned to R

  // ---- q_Nested == 0: A = I, matrix-free Woodbury (no n×n sparse) ----
  if(q_Nested == 0){
    arma::mat Ut_d = arma::mat(Ut);
    arma::mat U_d  = Ut_d.t();
    int Q = (int)Ut_d.n_rows;
    arma::mat IpUtU = arma::eye(Q, Q) + Ut_d * U_d;
    arma::mat M     = arma::inv_sympd(IpUtU);

    auto iV_vec = [&](const arma::vec& x) -> arma::vec {
      return x - U_d * (M * (Ut_d * x));
    };
    auto iV_mat_fn = [&](const arma::mat& A2) -> arma::mat {
      return A2 - U_d * (M * (Ut_d * A2));
    };

    arma::mat iVX  = iV_mat_fn(X);
    arma::mat denom = X.t() * iVX;                      // unscaled, for B
    B               = arma::solve(denom, X.t() * iV_vec(Y));
    H               = Y - X * B;
    HiVH            = arma::dot(H, iV_vec(H));
    s2resid         = REML ? HiVH / (n - p) : HiVH / n;

    // Form full n×n iV/s2resid for return: (I - U*M*Ut) / s2resid
    iV_scaled     = (arma::eye(n, n) - U_d * (M * Ut_d)) / s2resid;
    denom_scaled  = denom / s2resid;  // OPT-F: denom already computed O(n·p²)

  // ---- q_Nested > 0: Cholesky + triangular solves ----
  } else {
    arma::sp_mat A_sp = sp_mat(n, n); A_sp.eye();
    for(int j = 0; j < q_Nested; ++j){
      sp_mat nj = nested[j];
      A_sp = A_sp + pow(sn[j], 2) * nj;
    }
    arma::mat A1(A_sp);
    arma::mat chol_A;
    bool chol_ok = arma::chol(chol_A, A1);

    if(chol_ok){
      arma::mat Rlo = chol_A.t();   // lower L s.t. L*L' = A

      // Helper: solve A*out = RHS via stored Cholesky factors (O(n^2) per col)
      auto trisolve_m = [&](const arma::mat& RHS) -> arma::mat {
        arma::mat tmp = arma::solve(arma::trimatl(Rlo), RHS);
        return arma::solve(arma::trimatu(chol_A), tmp);
      };
      auto trisolve_v = [&](const arma::vec& rhs) -> arma::vec {
        arma::vec tmp = arma::solve(arma::trimatl(Rlo), rhs);
        return arma::solve(arma::trimatu(chol_A), tmp);
      };

      arma::mat iA_X = trisolve_m(X);
      arma::vec iA_Y = trisolve_v(Y);

      if(q_nonNested > 0){
        arma::mat U_d  = arma::mat(U);
        arma::mat Ut_d = arma::mat(Ut);
        int Q          = (int)Ut.n_rows;
        arma::mat iA_U = trisolve_m(U_d);
        arma::mat IpUtAU = arma::eye(Q, Q) + Ut_d * iA_U;
        arma::mat M      = arma::inv_sympd(IpUtAU);

        // iV*x = iA_x - iA_U*(M*(Ut*iA_x))
        arma::mat iVX   = iA_X - iA_U * (M * (Ut_d * iA_X));
        arma::vec iVY   = iA_Y - iA_U * (M * (Ut_d * iA_Y));
        arma::mat denom = X.t() * iVX;
        B               = arma::solve(denom, X.t() * iVY);
        H               = Y - X * B;
        arma::vec iA_H  = trisolve_v(H);
        arma::vec iVH   = iA_H - iA_U * (M * (Ut_d * iA_H));
        HiVH            = arma::dot(H, iVH);
        s2resid         = REML ? HiVH / (n - p) : HiVH / n;

        // OPT-A: form iA via chol2inv instead of n RHS triangular solves; then Woodbury for iV.
        // iA_U already computed above.
        arma::mat iA_full = phyr_chol2inv(chol_A);
        arma::mat iV_full = iA_full - iA_U * (M * iA_U.t());
        iV_scaled    = iV_full / s2resid;
        denom_scaled = denom / s2resid;  // OPT-F: skip O(n²·p) multiply

      } else {
        // q_nonNested == 0: iV = iA
        arma::mat denom = iA_X.t() * X;
        B               = arma::solve(denom, iA_X.t() * Y);
        H               = Y - X * B;
        arma::vec iA_H  = trisolve_v(H);
        HiVH            = arma::dot(H, iA_H);
        s2resid         = REML ? HiVH / (n - p) : HiVH / n;

        // OPT-A: chol2inv instead of n RHS solves; OPT-F: skip O(n²·p) multiply
        arma::mat iA_full = phyr_chol2inv(chol_A);
        iV_scaled    = iA_full / s2resid;
        denom_scaled = denom / s2resid;
      }

    } else {
      // Cholesky failed (A not PD for these params): LU fallback
      arma::mat iA1 = arma::inv(A1);
      arma::sp_mat iA = sp_mat(iA1);
      arma::sp_mat iV0;
      if(q_nonNested > 0){
        arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
        arma::mat IshUtAU = mat(Ishort + Ut * iA * U);
        iV0 = iA - iA * U * sp_mat(arma::inv_sympd(IshUtAU)) * Ut * iA;
      } else {
        iV0 = iA;
      }
      arma::mat iV(iV0);
      arma::mat denom = X.t() * iV * X;
      B       = arma::solve(denom, X.t() * iV * Y);
      H       = Y - X * B;
      HiVH    = arma::as_scalar(H.t() * iV * H);
      s2resid = REML ? HiVH / (n - p) : HiVH / n;
      iV_scaled    = iV / s2resid;
      denom_scaled = denom / s2resid;  // OPT-F
    }
  }

  rowvec s2r       = s2resid * pow(sr, 2);
  NumericVector s2n = s2resid * pow(sn, 2);
  arma::mat B_cov  = arma::inv_sympd(denom_scaled);
  arma::vec B_se   = arma::sqrt(B_cov.diag());

  return List::create(
    _["B"]       = B,
    _["B.se"]    = B_se,
    _["B.cov"]   = B_cov,
    _["sr"]      = sr,
    _["sn"]      = sn,
    _["s2n"]     = s2n,
    _["s2r"]     = s2r,
    _["s2resid"] = s2resid,
    _["iV"]      = iV_scaled,
    _["H"]       = H
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
# pglmm_gaussian_internal_cpp(par = s, X, Y, Zt = as(matrix(0, 0, 0), "sparseMatrix"),
#                             St = as(matrix(0, 0, 0), "sparseMatrix"), nested, REML,
#                             verbose, optimizer, maxit, 
#                             reltol, q, n, p, pi)
# res = pglmm_gaussian_internal_cpp(par = s, X, Y, Zt, St, nested, REML, 
#                             verbose, optimizer, maxit, 
#                             reltol, q, n, p, pi)
*/