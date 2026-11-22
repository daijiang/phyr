// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
#include "RcppArmadillo.h"
#include "pglmm_utils.h"   // detect_n_sp, verify_block_structure, block_chol
#include <nloptrAPI.h>     // OPT-D: NLopt C API via R_GetCCallable
#include <unordered_map>
#include <string>
#include <vector>

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
    n_sp_g = find_block_size_all(nested, q_Nested, n);
    n_site_g = (n_sp_g > 0) ? n / n_sp_g : 0;
    use_blk_g = (n_sp_g > 0);
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

    arma::mat IpUtU = arma::symmatu(arma::eye(Q, Q) + Ut_d * U_d);   // Q × Q
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
          arma::mat Ishort_Ut_iA_U = arma::symmatu(arma::eye(Q, Q) + Ut_iA_U);
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
          Ishort_Ut_iA_U_fb = arma::symmatu(mat(Ishort + Ut_iA_U));
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
    arma::mat IpUtU = arma::symmatu(arma::eye(Q, Q) + Ut_d * U_d);
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
        arma::mat IpUtAU = arma::symmatu(arma::eye(Q, Q) + Ut_d * iA_U);
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
        arma::mat IshUtAU = arma::symmatu(mat(Ishort + Ut * iA * U));
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
  arma::mat B_cov  = arma::inv_sympd(arma::symmatu(denom_scaled));
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
    _["H"]       = H
  );
}

// ---------------------------------------------------------------------------
// Models with non-nested random terms only (q_Nested == 0).
// Everything the profiled (RE)ML likelihood and its gradient need can be
// written with statistics that do not depend on the parameters:
//   G = Zt Zt' (Q x Q), ZX = Zt X (Q x p), Zy = Zt y (Q), X'X, X'y,
// which are computed once per fit. With c the Q-vector of standard deviations
// of the rows of Zt (Ut = diag(c) Zt, V = I + Ut' Ut), each evaluation needs
// only K = I + diag(c) G diag(c) (Q x Q): one Cholesky and one triangular
// inverse, with no dense work of size n (Woodbury identity and the matrix
// determinant lemma: V^-1 = I - Ut' K^-1 Ut, log|V| = log|K|).
struct GaussQ_pre {
  arma::mat G;      // Q x Q
  arma::mat ZX;     // Q x p
  arma::vec Zy;     // Q
  arma::mat XtX;    // p x p
  arma::vec Xty;    // p
  const arma::mat* X;
  const arma::vec* Y;
  arma::mat St_d;   // q x Q (0/1)
  int n, p, Q, q;
  bool REML;
};

static void gauss_q_precompute(GaussQ_pre& P, const arma::mat& X, const arma::vec& Y,
                               const arma::sp_mat& Zt, const arma::sp_mat& St, bool REML) {
  P.n = X.n_rows; P.p = X.n_cols; P.Q = Zt.n_rows; P.q = St.n_rows; P.REML = REML;
  // G = Zt Zt' from the unique columns of Zt (observations of the same levels
  // have the same column): G = Zu diag(count) Zu'
  {
    std::unordered_map<std::string, arma::uword> seen;
    std::vector<arma::uword> first;
    std::vector<double> cnt;
    for (int j = 0; j < P.n; ++j) {
      std::string key;
      for (arma::sp_mat::const_col_iterator it = Zt.begin_col(j); it != Zt.end_col(j); ++it) {
        arma::uword r = it.row();
        double v = *it;
        key.append((const char*)&r, sizeof(r));
        key.append((const char*)&v, sizeof(v));
      }
      std::unordered_map<std::string, arma::uword>::iterator f = seen.find(key);
      if (f == seen.end()) {
        seen.emplace(key, (arma::uword)first.size());
        first.push_back((arma::uword)j);
        cnt.push_back(1.0);
      } else {
        cnt[f->second] += 1.0;
      }
    }
    int n_u = (int)first.size();
    arma::mat Zu(P.Q, n_u, arma::fill::zeros);
    for (int g = 0; g < n_u; ++g) {
      for (arma::sp_mat::const_col_iterator it = Zt.begin_col(first[g]); it != Zt.end_col(first[g]); ++it) {
        Zu(it.row(), g) = *it;
      }
    }
    arma::mat ZuC = Zu;
    ZuC.each_row() %= arma::rowvec(cnt);
    P.G = ZuC * Zu.t();
  }
  P.ZX  = arma::mat(Zt * X);
  P.Zy  = arma::vec(Zt * Y);
  P.XtX = X.t() * X;
  P.Xty = X.t() * Y;
  P.X = &X;
  P.Y = &Y;
  P.St_d = arma::mat(St);
}

// Returns the same value as pglmm_gaussian_LL_cpp; fills grad if non-null, and
// optionally returns the pieces needed after the optimization.
static double gauss_q_ll_grad(const double* x, double* grad, const GaussQ_pre& P,
                              arma::vec* B_out = nullptr, double* HiVH_out = nullptr,
                              arma::mat* denom_out = nullptr) {
  int Q = P.Q, p = P.p, n = P.n, q = P.q;
  arma::rowvec sr(q);
  for (int r = 0; r < q; ++r) sr[r] = x[r];
  arma::vec c = arma::vectorise(sr * P.St_d);           // Q-vector

  // K = I + diag(c) G diag(c)
  arma::mat K = P.G;
  K.each_col() %= c;
  K.each_row() %= c.t();
  K.diag() += 1.0;
  arma::mat R;
  if (!arma::chol(R, arma::symmatu(K))) return 1e30;
  double logdetV = 2.0 * arma::accu(arma::log(R.diag()));
  arma::mat Rinv = arma::inv(arma::trimatu(R));          // K^-1 = Rinv * Rinv'

  arma::mat cZX = P.ZX; cZX.each_col() %= c;             // Ut X   (Q x p)
  arma::vec cZy = c % P.Zy;                              // Ut y   (Q)
  arma::mat TX  = Rinv * (Rinv.t() * cZX);               // K^-1 Ut X
  arma::vec Ty  = Rinv * (Rinv.t() * cZy);               // K^-1 Ut y
  arma::mat denom = P.XtX - cZX.t() * TX;                // X' iV X
  denom = arma::symmatu(denom);
  arma::vec num   = P.Xty - cZX.t() * Ty;                // X' iV y
  arma::vec B;
  if (!arma::solve(B, denom, num, arma::solve_opts::no_approx)) return 1e30;

  arma::vec H   = *P.Y - (*P.X) * B;                     // O(n p)
  arma::vec UtH = cZy - cZX * B;                         // Ut H  (Q)
  arma::vec b   = Ty - TX * B;                           // K^-1 Ut H
  double HiVH   = arma::dot(H, H) - arma::dot(UtH, b);

  int fac = P.REML ? (n - p) : n;
  double LL = 0.5 * (double(fac) * std::log(HiVH / double(fac)) + logdetV + double(fac));
  if (P.REML) {
    double ld, sg;
    arma::log_det(ld, sg, denom);
    LL += 0.5 * ld;
  }
  if (!std::isfinite(LL)) return 1e30;

  if (B_out) *B_out = B;
  if (HiVH_out) *HiVH_out = HiVH;
  if (denom_out) *denom_out = denom;

  if (grad) {
    arma::vec Mdiag = arma::sum(arma::square(Rinv), 1);  // diag(K^-1)
    arma::vec tc = (1.0 - Mdiag) - (double(fac) / HiVH) * (b % b);
    if (P.REML) {
      arma::mat UtF = arma::solve(denom, TX.t()).t();    // TX * denom^-1  (Q x p)
      tc -= arma::sum(TX % UtF, 1);
    }
    for (int k = 0; k < Q; ++k) tc[k] = (std::abs(c[k]) > 1e-15) ? tc[k] / c[k] : 0.0;
    for (int r = 0; r < q; ++r) {
      double g = arma::dot(P.St_d.row(r), tc.t());
      grad[r] = std::isfinite(g) ? g : 0.0;
    }
  }
  return LL;
}

// Estimates at the optimum, for models with non-nested random terms only
static List gauss_q_calc(const arma::vec& par, const GaussQ_pre& P) {
  int n = P.n, p = P.p, q = P.q;
  arma::vec B;
  double HiVH = 0.0;
  arma::mat denom;
  double LL = gauss_q_ll_grad(par.memptr(), nullptr, P, &B, &HiVH, &denom);
  if (LL >= 1e30) Rcpp::stop("Estimation of B failed.");
  double s2resid = P.REML ? HiVH / (n - p) : HiVH / n;
  arma::mat B_cov = arma::inv_sympd(arma::symmatu(denom / s2resid));
  arma::vec H = *P.Y - (*P.X) * B;
  arma::rowvec sr(q);
  for (int r = 0; r < q; ++r) sr[r] = par[r];
  arma::rowvec s2r = s2resid * arma::square(sr);
  return List::create(
    _["B"]       = B,
    _["B.se"]    = arma::vec(arma::sqrt(B_cov.diag())),
    _["B.cov"]   = B_cov,
    _["sr"]      = sr,
    _["sn"]      = NumericVector(0),
    _["s2n"]     = NumericVector(0),
    _["s2r"]     = s2r,
    _["s2resid"] = s2resid,
    _["H"]       = H
  );
}

// OPT-D: NLopt C API — parameter bundle and C-compatible objective.
// The callback calls pglmm_gaussian_LL_cpp directly (same TU, no R dispatch),
// eliminating the ~10-15ms per-call R→C++ roundtrip overhead of the nloptr R package.
struct GaussLL_ctx {
  const arma::mat*    X;
  const arma::vec*    Y;
  const arma::sp_mat* Zt;
  const arma::sp_mat* St;
  const Rcpp::List*   nested;
  bool REML;
  bool verbose;
  const GaussQ_pre*   pre;   // non-NULL for models with non-nested terms only
};

// OPT-J: combined LL + analytic gradient for L-BFGS.
//
// Key identity used throughout: for any Woodbury form
//   iV = iA - iA_U * M * (iA_U)',  M = (I + Ut*iA*U)^{-1},
// with Ut = iC*Zt and iC = diag(iC1), the diagonal of Zt*iV*Zt' satisfies
//   diag(Zt*iV*Zt')_k = (1 - M_kk) / iC1_k^2.
// This holds for all computation paths (non-nested, mixed, and—trivially—
// for nested-only paths where iV = iA and the M term is absent).
//
// Gradient of the concentrated LL w.r.t. theta_k:
//   dLL/d(sr_r) = dot(St_r, [(1-M_diag) - fac/HiVH*(b^2) - REML_UtF_norm2] / iC1)
//   dLL/d(sn_j) = sn_j * [tr(iV*nj) - fac/HiVH * v'*nj*v - REML tr(F'*nj*F)]
// where v = iVH, b = Ut*v, F = iVX*denom^{-1} (REML only),
//       fac = n-p (REML) or n (ML).
static double gauss_ll_and_grad(unsigned n_par, const double* x, double* grad,
                                const GaussLL_ctx& d) {
  const arma::mat&    X      = *d.X;
  const arma::vec&    Y      = *d.Y;
  const arma::sp_mat& Zt     = *d.Zt;
  const arma::sp_mat& St     = *d.St;
  const Rcpp::List&   nested = *d.nested;
  bool REML    = d.REML;

  int n = X.n_rows;
  int p = X.n_cols;
  int q_nonNested = (int)St.n_rows;
  int q_Nested    = (int)nested.size();
  int q_total     = q_nonNested + q_Nested;

  // ---- initialise gradient to 0 ----
  for (int k = 0; k < q_total; ++k) grad[k] = 0.0;

  // ---- extract parameters ----
  arma::vec iC1;       // Q-vector: iC1 = sr * St  (non-nested scaling)
  arma::mat Ut_d;      // Q × n dense
  arma::mat U_d;       // n × Q dense
  arma::mat St_d;      // q_nonNested × Q dense (for gradient assembly)
  int Q = 0;

  if (q_nonNested > 0) {
    arma::rowvec sr(q_nonNested);
    for (int r = 0; r < q_nonNested; ++r) sr[r] = x[r];
    arma::mat iC0 = sr * arma::mat(St);
    iC1  = arma::vectorise(iC0, 0);      // Q-vector
    Q    = (int)iC1.n_elem;
    arma::sp_mat iC_sp = arma::sp_mat(arma::diagmat(iC1));
    arma::sp_mat Ut_sp = iC_sp * Zt;
    Ut_d = arma::mat(Ut_sp);             // Q × n
    U_d  = Ut_d.t();                     // n × Q
    St_d = arma::mat(St);                // q_nonNested × Q
  }

  arma::vec sn_vec(q_Nested);
  for (int j = 0; j < q_Nested; ++j) sn_vec[j] = x[q_nonNested + j];

  // block-diagonal detection (same logic as LL function)
  int  n_sp_g = 0, n_site_g = 0;
  bool use_blk_g = false;
  if (q_Nested > 0) {
    n_sp_g = find_block_size_all(nested, q_Nested, n);
    n_site_g = (n_sp_g > 0) ? n / n_sp_g : 0;
    use_blk_g = (n_sp_g > 0);
  }

  // ---- branch-local state (filled per path) ----
  double logdetV = 0.0, HiVH = 0.0;
  arma::mat denom;     // p × p
  arma::vec v;         // n-vector: iVH (= iV * H, key gradient intermediate)

  // Non-nested Woodbury M (for grad via 1-M_diag):
  arma::mat M_mat;          // Q × Q, set when Woodbury is used
  bool have_M = false;

  // For nested trace terms:
  // Path B: block Cholesky factors
  std::vector<arma::mat> chols_blk_g;
  // Path C: full Cholesky factor
  arma::mat chol_A_g;
  bool have_chol_A = false;
  // Path C2 (mixed): iA*U and Ut*iA*U
  arma::mat iA_U_g;    // n × Q
  arma::mat ZtiAU_g;   // Q × Q = Ut*iA*U

  // For REML correction: iVX (n × p), used to build F = iVX * denom^{-1}
  arma::mat iVX_d;     // set in each branch
  // LU fallback: retain the materialised dense iV for nested trace computation
  arma::mat iV_lu;     // only set on LU fallback path

  // ===== PATH A: non-nested only, matrix-free Woodbury =====
  if (q_Nested == 0) {
    arma::mat IpUtU = arma::symmatu(arma::eye(Q, Q) + Ut_d * U_d);
    arma::mat M     = arma::inv_sympd(IpUtU);
    double signV;
    arma::log_det(logdetV, signV, IpUtU);
    if (!std::isfinite(logdetV)) {
      arma::mat lgm = arma::chol(IpUtU);
      logdetV = 2.0 * arma::accu(arma::log(lgm.diag()));
    }

    auto iV_vec = [&](const arma::vec& z) -> arma::vec {
      return z - U_d * (M * (Ut_d * z));
    };
    auto iV_mat_fn = [&](const arma::mat& A2) -> arma::mat {
      return A2 - U_d * (M * (Ut_d * A2));
    };

    iVX_d        = iV_mat_fn(X);
    denom        = X.t() * iVX_d;
    arma::vec B  = arma::solve(denom, X.t() * iV_vec(Y));
    arma::vec H  = Y - X * B;
    v            = iV_vec(H);
    HiVH         = arma::dot(H, v);

    M_mat   = M;
    have_M  = true;

  // ===== PATH B: nested-only, block Cholesky =====
  } else if (q_nonNested == 0 && use_blk_g) {
    std::vector<arma::mat> nj_fixed(n_site_g,
                                    arma::mat(n_sp_g, n_sp_g, arma::fill::zeros));
    for (int j = 0; j < q_Nested; ++j) {
      double snj2 = sn_vec[j] * sn_vec[j];
      arma::sp_mat njsp = nested[j];
      for (int k = 0; k < n_site_g; ++k) {
        int s = k * n_sp_g;
        nj_fixed[k] += snj2 *
          arma::mat(njsp(arma::span(s, s+n_sp_g-1), arma::span(s, s+n_sp_g-1)));
      }
    }
    chols_blk_g.resize(n_site_g);
    arma::vec pq_ones(n, arma::fill::ones);
    double logdetA = block_chol(chols_blk_g, pq_ones, nj_fixed, n_sp_g, n_site_g);

    auto blk_solve_v = [&](const arma::vec& rhs) -> arma::vec {
      arma::vec out(n);
      for (int k = 0; k < n_site_g; ++k) {
        int s = k * n_sp_g;
        arma::mat Rlo_k = chols_blk_g[k].t();
        arma::vec tmp = arma::solve(arma::trimatl(Rlo_k), rhs.subvec(s, s+n_sp_g-1));
        out.subvec(s, s+n_sp_g-1) = arma::solve(arma::trimatu(chols_blk_g[k]), tmp);
      }
      return out;
    };
    auto blk_solve_m = [&](const arma::mat& RHS) -> arma::mat {
      arma::mat out(RHS.n_rows, RHS.n_cols);
      for (int k = 0; k < n_site_g; ++k) {
        int s = k * n_sp_g;
        arma::mat Rlo_k = chols_blk_g[k].t();
        arma::mat tmp = arma::solve(arma::trimatl(Rlo_k), RHS.rows(s, s+n_sp_g-1));
        out.rows(s, s+n_sp_g-1) = arma::solve(arma::trimatu(chols_blk_g[k]), tmp);
      }
      return out;
    };

    arma::mat iA_X = blk_solve_m(X);
    iVX_d          = iA_X;
    denom          = iA_X.t() * X;
    arma::vec B    = arma::solve(denom, iA_X.t() * Y);
    arma::vec H    = Y - X * B;
    v              = blk_solve_v(H);
    HiVH           = arma::dot(H, v);
    logdetV        = logdetA;

  // ===== PATH B2: nested (block-diagonal) + non-nested, block Cholesky + Woodbury =====
  } else if (use_blk_g) {
    std::vector<arma::mat> nj_fixed(n_site_g,
                                    arma::mat(n_sp_g, n_sp_g, arma::fill::zeros));
    for (int j = 0; j < q_Nested; ++j) {
      double snj2 = sn_vec[j] * sn_vec[j];
      arma::sp_mat njsp = nested[j];
      for (int k = 0; k < n_site_g; ++k) {
        int s = k * n_sp_g;
        nj_fixed[k] += snj2 *
          arma::mat(njsp(arma::span(s, s+n_sp_g-1), arma::span(s, s+n_sp_g-1)));
      }
    }
    chols_blk_g.resize(n_site_g);
    arma::vec pq_ones(n, arma::fill::ones);
    double logdetA = block_chol(chols_blk_g, pq_ones, nj_fixed, n_sp_g, n_site_g);
    auto blk_solve_m = [&](const arma::mat& RHS) -> arma::mat {
      arma::mat out(RHS.n_rows, RHS.n_cols);
      for (int k = 0; k < n_site_g; ++k) {
        int s = k * n_sp_g;
        arma::mat tmp = arma::solve(arma::trimatl(chols_blk_g[k].t()), RHS.rows(s, s+n_sp_g-1));
        out.rows(s, s+n_sp_g-1) = arma::solve(arma::trimatu(chols_blk_g[k]), tmp);
      }
      return out;
    };
    arma::mat iA_X = blk_solve_m(X);
    arma::mat iA_U = blk_solve_m(U_d);
    iA_U_g = iA_U;
    arma::mat ZtiAU = Ut_d * iA_U;
    ZtiAU_g = ZtiAU;
    int Qc = (int)Ut_d.n_rows;
    arma::mat Ishort = arma::symmatu(arma::eye(Qc, Qc) + ZtiAU);
    arma::mat M      = arma::inv_sympd(Ishort);
    M_mat = M; have_M = true;
    arma::mat Ut_iA_X = Ut_d * iA_X;
    arma::mat iA_Y    = blk_solve_m(arma::mat(Y));
    arma::vec Ut_iA_Y = Ut_d * iA_Y.col(0);
    denom            = iA_X.t() * X - Ut_iA_X.t() * M * Ut_iA_X;
    arma::vec B      = arma::solve(denom, iA_X.t() * Y - Ut_iA_X.t() * M * Ut_iA_Y);
    arma::vec H      = Y - X * B;
    arma::mat iA_H   = blk_solve_m(arma::mat(H));
    arma::vec Ut_iA_H = Ut_d * iA_H.col(0);
    arma::vec M_Ut_iA_H = M * Ut_iA_H;
    v                = iA_H.col(0) - iA_U * M_Ut_iA_H;
    HiVH             = arma::dot(H, iA_H.col(0)) - arma::as_scalar(Ut_iA_H.t() * M_Ut_iA_H);
    iVX_d = iA_X - iA_U * (M * Ut_iA_X);
    double signV;
    arma::log_det(logdetV, signV, Ishort);
    logdetV += logdetA;

  // ===== PATH C: full dense path =====
  } else {
    arma::sp_mat A_sp = arma::sp_mat(n, n); A_sp.eye();
    for (int j = 0; j < q_Nested; ++j)
      A_sp = A_sp + sn_vec[j]*sn_vec[j] * arma::sp_mat(nested[j]);
    arma::mat A1(A_sp);

    arma::mat chol_A;
    bool chol_ok = arma::chol(chol_A, A1);

    if (chol_ok) {
      double logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
      arma::mat Rlo  = chol_A.t();
      chol_A_g    = chol_A;
      have_chol_A = true;

      arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
      iA_X           = arma::solve(arma::trimatu(chol_A), iA_X);

      if (q_nonNested > 0) {
        // Woodbury on top of Cholesky
        arma::mat iA_U = arma::solve(arma::trimatl(Rlo), U_d);
        iA_U           = arma::solve(arma::trimatu(chol_A), iA_U);
        iA_U_g         = iA_U;  // retain for gradient

        arma::mat ZtiAU = Ut_d * iA_U;          // Q × Q
        ZtiAU_g         = ZtiAU;
        int Qc = (int)Ut_d.n_rows;
        arma::mat Ishort = arma::symmatu(arma::eye(Qc, Qc) + ZtiAU);
        arma::mat M      = arma::inv_sympd(Ishort);
        M_mat   = M;
        have_M  = true;

        arma::mat Ut_iA_X = Ut_d * iA_X;
        arma::mat iA_Y    = arma::solve(arma::trimatl(Rlo), arma::mat(Y));
        iA_Y              = arma::solve(arma::trimatu(chol_A), iA_Y);
        arma::vec Ut_iA_Y = Ut_d * iA_Y.col(0);

        denom            = iA_X.t() * X - Ut_iA_X.t() * M * Ut_iA_X;
        arma::vec B      = arma::solve(denom, iA_X.t() * Y - Ut_iA_X.t() * M * Ut_iA_Y);
        arma::vec H      = Y - X * B;

        arma::mat iA_H   = arma::solve(arma::trimatl(Rlo), arma::mat(H));
        iA_H             = arma::solve(arma::trimatu(chol_A), iA_H);
        arma::vec Ut_iA_H = Ut_d * iA_H.col(0);
        arma::vec M_Ut_iA_H = M * Ut_iA_H;
        v                = iA_H.col(0) - iA_U * M_Ut_iA_H;
        HiVH             = arma::dot(H, iA_H.col(0)) -
                           arma::as_scalar(Ut_iA_H.t() * M_Ut_iA_H);

        // iVX for REML = iA_X - iA_U * M * Ut_iA_X
        iVX_d = iA_X - iA_U * (M * Ut_iA_X);

        double signV;
        arma::log_det(logdetV, signV, Ishort);
        logdetV += logdetA;

      } else {
        // nested-only, full dense Cholesky
        iVX_d          = iA_X;
        denom          = iA_X.t() * X;
        arma::vec B    = arma::solve(denom, iA_X.t() * Y);
        arma::vec H    = Y - X * B;
        arma::mat iA_H = arma::solve(arma::trimatl(Rlo), arma::mat(H));
        iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
        v              = iA_H.col(0);
        HiVH           = arma::dot(H, v);
        logdetV        = logdetA;
      }

    } else {
      // LU fallback: materialise dense iV (all-dense arithmetic)
      arma::mat iA1 = arma::inv(A1);
      double logdetA_fb, sgn;
      arma::log_det(logdetA_fb, sgn, A1);

      arma::mat iV_dense;
      if (q_nonNested > 0) {
        arma::mat iA_U_fb = iA1 * U_d;
        arma::mat ZtiAU   = Ut_d * iA_U_fb;
        int Qlu = (int)Ut_d.n_rows;
        arma::mat Ishort  = arma::symmatu(arma::eye(Qlu, Qlu) + ZtiAU);
        arma::mat M       = arma::inv_sympd(Ishort);
        M_mat   = M;
        have_M  = true;
        iA_U_g  = iA_U_fb;
        ZtiAU_g = ZtiAU;
        iV_dense = iA1 - iA_U_fb * M * iA_U_fb.t();
        double signV;
        arma::log_det(logdetV, signV, Ishort);
        logdetV += logdetA_fb;
      } else {
        iV_dense = iA1;
        logdetV  = logdetA_fb;
      }
      iV_lu  = iV_dense;   // save for nested trace in gradient section
      denom  = X.t() * iV_dense * X;
      arma::vec B = arma::solve(denom, X.t() * iV_dense * Y);
      arma::vec H = Y - X * B;
      v      = iV_dense * H;
      HiVH   = arma::dot(H, v);
      iVX_d  = iV_dense * X;
    }
  }

  // ---- Compute LL ----
  double LL;
  int factor = REML ? (n - p) : n;
  if (REML) {
    double s2c = HiVH / double(factor);
    double ldL, sgnL;
    arma::log_det(ldL, sgnL, denom);
    LL = 0.5 * (double(factor) * std::log(s2c) + logdetV + double(factor) + ldL);
  } else {
    double s2c = HiVH / double(factor);
    LL = 0.5 * (double(factor) * std::log(s2c) + logdetV + double(factor));
  }
  if (!std::isfinite(LL)) return 1e30;

  // ---- Build REML correction: F = iVX * denom^{-1} ----
  arma::mat F_mat;
  arma::mat UtF;          // Q × p (for non-nested REML term)
  arma::vec rowNorm2_UtF; // Q-vector
  if (REML && q_nonNested > 0) {
    arma::mat denom_inv = arma::solve(denom, arma::eye(p, p));
    F_mat = iVX_d * denom_inv;
    UtF   = Ut_d * F_mat;
    // Correct quadratic form: u_k' * P * u_k = UtF[k,:] * denom * UtF[k,:]'
    // where P = iVX*(X'iVX)^{-1}*X'iV. Using ||UtF[k,:]||^2 (identity instead
    // of denom) is wrong whenever X'iVX != I_p.
    rowNorm2_UtF = arma::sum((UtF * denom) % UtF, 1);
  }
  arma::mat F_mat_nested;  // n × p, for nested REML term (computed on demand)
  bool F_nested_built = false;
  if (REML && q_Nested > 0) {
    if (q_nonNested > 0) {
      F_mat_nested = F_mat;  // already built above
    } else {
      arma::mat denom_inv = arma::solve(denom, arma::eye(p, p));
      F_mat_nested = iVX_d * denom_inv;
    }
    F_nested_built = true;
  }

  // ---- Gradient: non-nested parameters sr_r ----
  if (q_nonNested > 0 && have_M) {
    // diag(Zt*iV*Zt')_k * iC1_k^2 = (ZtiAU*M)_kk for general A=I+nested.
    // When A=I (Path A, q_Nested==0): ZtiAU*M = (Ut*U)*(I+Ut*U)^{-1} = I-M,
    // so the diagonal simplifies to 1-M_kk and ZtiAU_g is not stored.
    // When A!=I (Path C2 or LU fallback with mixed): ZtiAU_g = Ut*iA*U is stored.
    arma::vec trace_diag;
    if (ZtiAU_g.n_rows == (arma::uword)Q && Q > 0) {
      arma::mat ZtiAU_M = ZtiAU_g * M_mat;
      trace_diag = ZtiAU_M.diag();             // Path C2 / LU fallback mixed
    } else {
      trace_diag = 1.0 - M_mat.diag();         // Path A (A=I)
    }

    arma::vec b      = Ut_d * v;                  // Q-vector: Ut * iVH

    arma::vec tc_raw = trace_diag
                       - (double(factor) / HiVH) * (b % b)
                       - (REML ? rowNorm2_UtF : arma::zeros<arma::vec>(Q));

    // guard against near-zero iC1 (degenerate sr = 0)
    arma::vec tc_scaled(Q);
    for (int k = 0; k < Q; ++k)
      tc_scaled[k] = (std::abs(iC1[k]) > 1e-15) ? tc_raw[k] / iC1[k] : 0.0;

    for (int r = 0; r < q_nonNested; ++r)
      grad[r] = arma::dot(St_d.row(r), tc_scaled);
  }

  // ---- Gradient: nested parameters sn_j ----
  if (q_Nested > 0) {
    for (int j = 0; j < q_Nested; ++j) {
      arma::sp_mat nj_sp = nested[j];

      // tr(iV * nested[j]) — path-dependent
      double tr_iV_nj = 0.0;

      if (!chols_blk_g.empty()) {
        // PATH B / B2: block Cholesky; the Woodbury correction is done by block:
        //   tr(M iA_U' N iA_U) = sum_k accu((N_k iA_U[k,]) % (iA_U M)[k,])
        bool wood = (q_nonNested > 0 && have_M);
        arma::mat PM;
        if (wood) PM = iA_U_g * M_mat;                       // n x Q
        for (int k = 0; k < n_site_g; ++k) {
          int s = k * n_sp_g;
          arma::mat nj_blk = arma::mat(
            nj_sp(arma::span(s, s+n_sp_g-1), arma::span(s, s+n_sp_g-1)));
          arma::mat Rlo_k = chols_blk_g[k].t();
          arma::mat W     = arma::solve(arma::trimatl(Rlo_k), nj_blk);
          W               = arma::solve(arma::trimatu(chols_blk_g[k]), W);
          tr_iV_nj       += arma::trace(W);
          if (wood)
            tr_iV_nj -= arma::accu((nj_blk * iA_U_g.rows(s, s+n_sp_g-1)) % PM.rows(s, s+n_sp_g-1));
        }
      } else if (have_chol_A) {
        // PATH C: full dense Cholesky; iV = iA (nested-only) or Woodbury (mixed)
        // tr(iA * nested[j]) via triangular solve
        arma::mat Rlo  = chol_A_g.t();
        arma::mat Nj_d = arma::mat(nj_sp);
        arma::mat W    = arma::solve(arma::trimatl(Rlo), Nj_d);
        W              = arma::solve(arma::trimatu(chol_A_g), W);
        double tr_iA_nj = arma::trace(W);

        if (q_nonNested > 0 && have_M) {
          // Woodbury correction: -tr(M * (iA_U)' * nested[j] * iA_U)
          arma::mat NjiAU   = nj_sp * iA_U_g;       // n × Q  (sparse * dense)
          double tr_corr    = arma::accu((iA_U_g.t() * NjiAU) % M_mat);
          tr_iV_nj = tr_iA_nj - tr_corr;
        } else {
          tr_iV_nj = tr_iA_nj;
        }
      } else {
        // LU fallback: iV was materialised and saved in iV_lu.
        // tr(iV * nested[j]) = Frobenius inner product (both symmetric).
        tr_iV_nj = arma::accu(iV_lu % arma::mat(nj_sp));
      }

      // v' * nested[j] * v  (sparse matvec)
      double vtnjv = arma::dot(v, nj_sp * v);

      // REML: tr(P * nj) = tr(iVX' * nj * F) = accu(iVX_d % (nj * F))
      // where P = iVX*(X'iVX)^{-1}*X'iV and F = iVX*(X'iVX)^{-1}.
      // Using tr(F'*nj*F) is wrong: it has (X'iVX)^{-2} instead of (X'iVX)^{-1}.
      double reml_nj = 0.0;
      if (REML && F_nested_built)
        reml_nj = arma::accu(iVX_d % (nj_sp * F_mat_nested));

      grad[q_nonNested + j] =
        sn_vec[j] * (tr_iV_nj - double(factor) / HiVH * vtnjv - reml_nj);
    }
  }

  // check gradient finiteness
  for (int k = 0; k < q_total; ++k)
    if (!std::isfinite(grad[k])) { grad[k] = 0.0; }

  return LL;
}

static double gauss_ll_nlopt_cb(unsigned n_par, const double* x,
                                double* grad, void* f_data) {
  const GaussLL_ctx* d = static_cast<const GaussLL_ctx*>(f_data);
  if (d->pre != nullptr) {
    try {
      return gauss_q_ll_grad(x, grad, *d->pre);
    } catch (...) {
      if (grad != nullptr) for (unsigned k = 0; k < n_par; ++k) grad[k] = 0.0;
      return 1e30;
    }
  }
  if (grad != nullptr) {
    try {
      return gauss_ll_and_grad(n_par, x, grad, *d);
    } catch (...) {
      for (unsigned k = 0; k < n_par; ++k) grad[k] = 0.0;
      return 1e30;
    }
  }
  try {
    Rcpp::NumericVector par(x, x + (int)n_par);
    return pglmm_gaussian_LL_cpp(par, *d->X, *d->Y, *d->Zt, *d->St,
                                  *d->nested, d->REML, d->verbose);
  } catch (...) {
    return 1e30;  // non-PD or other failure: penalize so optimizer moves away
  }
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

  arma::vec par_opt0;
  double    LL;
  int       convcode;
  arma::vec niter;

  // statistics that do not depend on the parameters, computed once
  // (models with non-nested random terms only)
  bool use_q = (nested.size() == 0 && St.n_rows > 0 && Zt.n_rows > 0);
  GaussQ_pre pre_q;
  if (use_q) gauss_q_precompute(pre_q, X, Y, Zt, St, REML);

  if (optimizer == "Nelder-Mead") {
    // R optim path (Nelder-Mead / L-BFGS-B): unchanged from original.
    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim_fn = stats["optim"];
    Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
    Rcpp::Function pglmm_gaussian_LL_cpp_fxn = phyr_pkg["pglmm_gaussian_LL_cpp"];

    Rcpp::List opt;
    if (q > 1) {
      opt = optim_fn(_["par"] = par, _["fn"] = pglmm_gaussian_LL_cpp_fxn,
                     _["X"] = X, _["Y"] = Y, _["Zt"] = Zt, _["St"] = St,
                     _["nested"] = nested, _["REML"] = REML, _["verbose"] = verbose,
                     _["method"] = "Nelder-Mead",
                     _["control"] = List::create(_["maxit"] = maxit, _["reltol"] = reltol));
    } else {
      opt = optim_fn(_["par"] = par, _["fn"] = pglmm_gaussian_LL_cpp_fxn,
                     _["X"] = X, _["Y"] = Y, _["Zt"] = Zt, _["St"] = St,
                     _["nested"] = nested, _["REML"] = REML, _["verbose"] = verbose,
                     _["method"] = "L-BFGS-B",
                     _["control"] = List::create(_["maxit"] = maxit));
    }
    par_opt0 = abs(real(as<arma::vec>(opt["par"])));
    LL       = as_scalar(as<double>(opt["value"]));
    convcode = as<int>(opt["convergence"]);
    niter    = as<arma::vec>(opt["counts"]);

  } else {
    // OPT-D / OPT-J: NLopt C API.
    // Derivative-free: BOBYQA, Nelder-Mead, SBPLX.
    // Gradient-based:  LBFGS (OPT-J: analytic gradient, ~5-10x fewer iterations).
    nlopt_algorithm algor;
    if      (optimizer == "bobyqa")            algor = NLOPT_LN_BOBYQA;
    else if (optimizer == "nelder-mead-nlopt") algor = NLOPT_LN_NELDERMEAD;
    else if (optimizer == "lbfgs")             algor = NLOPT_LD_LBFGS;
    else                                       algor = NLOPT_LN_SBPLX;

    GaussLL_ctx ctx = {&X, &Y, &Zt, &St, &nested, REML, verbose, use_q ? &pre_q : nullptr};
    nlopt_opt opt_c = nlopt_create(algor, (unsigned)q);
    nlopt_set_min_objective(opt_c, gauss_ll_nlopt_cb, &ctx);
    nlopt_set_ftol_rel(opt_c, reltol);
    nlopt_set_ftol_abs(opt_c, reltol);
    nlopt_set_xtol_rel(opt_c, 0.0001);
    nlopt_set_maxeval(opt_c, maxit);

    // No lower bounds needed for Gaussian: the concentrated LL is symmetric in
    // the sign of each sd parameter (V depends on sd², so LL(x) = LL(-x)).
    // The optimizer can roam freely in ℝ^q; abs() applied to the solution below
    // always recovers the correct non-negative sd estimate.
    // (Contrast: binary pql_ll_and_grad uses abs(x) internally but returns the
    //  unsigned gradient, causing sign errors for x < 0 → bounds needed there.)

    std::vector<double> x_vec(par.begin(), par.end());
    double opt_val = 1e30;
    nlopt_result res = nlopt_optimize(opt_c, x_vec.data(), &opt_val);
    convcode = (int)res;
    nlopt_destroy(opt_c);

    par_opt0 = arma::abs(arma::vec(x_vec.data(), (arma::uword)q));
    LL       = opt_val;
    // nlopt_get_numevals not exposed via nloptrAPI.h; return placeholder.
    niter    = arma::vec(1, arma::fill::ones);
  }

  NumericVector par_opt = wrap(par_opt0);

  // calculate coef
  List out = use_q ? gauss_q_calc(par_opt0, pre_q)
                   : pglmm_gaussian_LL_calc_cpp(par_opt, X, Y, Zt, St, nested, REML);
  // restricted log-likelihood as in lme4, glmmTMB and nlme (the likelihood
  // integrated over the fixed effects, without the term 0.5 log|X'X| of the
  // density of orthonormal error contrasts)
  double logLik;
  if (REML) {
    logLik = -0.5 * (n - p) * log(2 * Pi) - LL;
  } else {
    logLik = -0.5 * n * log(2 * Pi) - LL;
  }

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