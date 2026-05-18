// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
#include "RcppArmadillo.h"
#include "pglmm_utils.h"   // detect_n_sp, verify_block_structure, block_chol, phyr_chol2inv
#include <nloptrAPI.h>     // OPT-K: NLopt C API for PQL analytic gradient
using namespace Rcpp;
using namespace arma;

// ============================================================
// Internal helpers
// ============================================================

// Build Ut (Q×n sparse) from non-nested variance params + design matrices
static arma::sp_mat build_Ut_helper(const NumericVector& par,
                                    const arma::sp_mat& St,
                                    const arma::sp_mat& Zt,
                                    int q_nn) {
  IntegerVector idx = seq_len(q_nn) - 1;
  NumericVector sr0 = par[idx];
  rowvec sr = real(as<rowvec>(sr0));
  arma::vec  iC1 = vectorise(sr * St, 0);
  arma::sp_mat iC = sp_mat(diagmat(iC1));
  return iC * Zt;
}

// ============================================================
// pglmm_iV_logdetV_cpp  (unchanged public API)
// ============================================================
// [[Rcpp::export]]
List pglmm_iV_logdetV_cpp(NumericVector par, arma::vec mu,
                           const arma::sp_mat& Zt, const arma::sp_mat& St,
                           const List& nested, bool logdet,
                           const std::string family, arma::vec totalSize){
  int q_nn = St.n_rows;
  arma::sp_mat Ut, U;
  if(q_nn > 0){
    IntegerVector idx = seq_len(q_nn) - 1;
    NumericVector sr0 = par[idx];
    rowvec sr = real(as<rowvec>(sr0));
    arma::vec iC1 = vectorise(sr * St, 0);
    Ut = sp_mat(diagmat(iC1)) * Zt;
    U  = trans(Ut);
  }
  int q_Nested = nested.size();
  NumericVector sn;
  if(q_Nested > 0){
    IntegerVector idx2 = wrap(seq(q_nn, q_nn + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    rowvec sn1 = real(as<rowvec>(sn0));
    sn = as<NumericVector>(wrap(sn1));
  }

  arma::sp_mat iV0;
  arma::mat Ishort_Ut_iA_U;
  double logdetV = 0.0, signV, logdetiA, signiA;

  if(q_Nested == 0){
    arma::vec pq = totalSize % mu % (1 - mu);
    arma::sp_mat iA;
    if(family == "binomial") iA = sp_mat(diagmat(pq));
    if(family == "poisson")  iA = sp_mat(diagmat(mu));
    arma::sp_mat Ishort(Ut.n_rows, Ut.n_rows); Ishort.eye();
    Ishort_Ut_iA_U = mat(Ishort + Ut * iA * U);
    arma::mat i_Ish = inv(Ishort_Ut_iA_U);
    iV0 = iA - iA * U * sp_mat(i_Ish) * Ut * iA;
    if(logdet){
      log_det(logdetV, signV, Ishort_Ut_iA_U);
      log_det(logdetiA, signiA, mat(iA));
      logdetV -= logdetiA;
      if(std::isinf(logdetV)){
        arma::mat chol_Ish = arma::chol(Ishort_Ut_iA_U);
        logdetV = 2*sum(log(chol_Ish.diag())) - logdetiA;
      }
    }
  } else {
    arma::vec pq = 1.0 / (totalSize % mu % (1 - mu));
    arma::sp_mat A;
    if(family == "binomial") A = sp_mat(diagmat(pq));
    if(family == "poisson")  A = sp_mat(diagmat(1.0/mu));
    for(int j = 0; j < q_Nested; ++j){
      double snj2 = pow((double)sn[j], 2);
      sp_mat nj = nested[j];   // implicit conversion — do not use as<>
      A = A + snj2 * nj;
    }
    arma::mat A1(A);
    arma::mat chol_A; double logdetA = 0.0;
    arma::sp_mat iA;
    if(arma::chol(chol_A, A1)){
      logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
      // OPT-E: chol2inv instead of inv(trimatu) + multiply
      iA = sp_mat(phyr_chol2inv(chol_A));
    } else {
      double sgn; arma::log_det(logdetA, sgn, A1);
      iA = sp_mat(arma::inv(A1));
    }
    if(q_nn > 0){
      arma::sp_mat Ishort(Ut.n_rows, Ut.n_rows); Ishort.eye();
      Ishort_Ut_iA_U = mat(Ishort + Ut * iA * U);
      iV0 = iA - iA * U * sp_mat(arma::inv_sympd(Ishort_Ut_iA_U)) * Ut * iA;
    } else {
      iV0 = iA;
    }
    if(logdet){
      if(q_nn > 0){
        double s2; log_det(logdetV, s2, Ishort_Ut_iA_U);
        logdetV += logdetA;
      } else {
        logdetV = logdetA;
      }
    }
  }
  if(logdet) return List::create(_["iV"]=iV0, _["logdetV"]=logdetV);
  return List::create(_["iV"]=iV0);
}

// [[Rcpp::export]]
arma::sp_mat pglmm_V(NumericVector par, const arma::sp_mat& Zt,
                     const arma::sp_mat& St, arma::vec mu,
                     const List& nested, bool missing_mu,
                     const std::string family, arma::vec totalSize){
  int q_nn = St.n_rows;
  arma::sp_mat Ut, U;
  if(q_nn > 0){
    IntegerVector idx = seq_len(q_nn) - 1;
    NumericVector sr0 = par[idx];
    rowvec sr = real(as<rowvec>(sr0));
    arma::vec iC1 = vectorise(sr * St, 0);
    Ut = sp_mat(diagmat(iC1)) * Zt;
    U  = trans(Ut);
  }
  int q_Nested = nested.size();
  NumericVector sn;
  if(q_Nested > 0){
    IntegerVector idx2 = wrap(seq(q_nn, q_nn + q_Nested - 1));
    NumericVector sn0 = par[idx2];
    rowvec sn1 = real(as<rowvec>(sn0));
    sn = as<NumericVector>(wrap(sn1));
  }
  arma::mat iW;
  if(missing_mu){
    iW = mat(Zt.n_cols, Zt.n_cols, fill::zeros);
  } else {
    arma::vec pq = 1.0/(totalSize % mu % (1 - mu));
    if(family == "binomial") iW = diagmat(pq);
    if(family == "poisson")  iW = diagmat(1.0/mu);
  }
  arma::sp_mat A = sp_mat(iW);
  for(int j = 0; j < q_Nested; ++j){
    double snj2 = pow((double)sn[j], 2);
    sp_mat nj = nested[j];
    A = A + snj2 * nj;
  }
  return (q_nn > 0) ? A + U * Ut : A;
}

// ============================================================
// pglmm_LL_cpp: OPT5a (q_Nested==0) + OPT6a (block-diagonal A)
// ============================================================
// [[Rcpp::export]]
double pglmm_LL_cpp(NumericVector par, const arma::vec& H,
                    const arma::mat& X, const arma::sp_mat& Zt,
                    const arma::sp_mat& St, const arma::vec& mu,
                    const List& nested, bool REML, bool verbose,
                    const std::string family, arma::vec totalSize){
  par = abs(par);
  int q_nn     = St.n_rows;
  int q_Nested = (int)nested.size();
  double LL    = 0.0;

  // ---- OPT5a: q_Nested==0, diagonal iA ----
  if(q_Nested == 0 && q_nn > 0){
    arma::sp_mat Ut_sp = build_Ut_helper(par, St, Zt, q_nn);
    arma::mat Ut_d = arma::mat(Ut_sp);
    arma::mat U_d  = Ut_d.t();
    int Q = (int)Ut_d.n_rows;

    arma::vec W;
    if(family=="binomial") W = totalSize % mu % (1.0-mu);
    if(family=="poisson")  W = mu;

    arma::mat WU = U_d; WU.each_col() %= W;
    arma::mat IpUtWU = arma::eye(Q,Q) + Ut_d * WU;
    arma::mat M = arma::inv_sympd(IpUtWU);

    double lsign, logdetV;
    arma::log_det(logdetV, lsign, IpUtWU);
    logdetV -= arma::accu(arma::log(W));

    arma::vec WH  = W % H;
    arma::vec iVH = WH - WU * (M * (Ut_d * WH));

    if(REML){
      arma::mat WX = X; WX.each_col() %= W;
      arma::mat iVX = WX - WU * (M * (Ut_d * WX));
      double ld, sg; arma::log_det(ld, sg, X.t() * iVX);
      LL = 0.5*(logdetV + arma::dot(H, iVH) + ld);
    } else {
      LL = 0.5*(logdetV + arma::dot(H, iVH));
    }

  // ---- q_Nested>0: Gaussian-style triangular-solve path (no inv, no sp_mat, no List) ----
  } else if(q_Nested > 0){
    // Build A = diag(pq) + sum_j snj^2 * nj exactly as pglmm_iV_logdetV_cpp does,
    // so logdetV is bit-identical to the original; then use triangular solves instead
    // of inv()+sp_mat() to eliminate the redundant O(n^3) inversion.
    arma::vec pq;
    if(family=="binomial") pq = 1.0/(totalSize % mu % (1.0-mu));
    if(family=="poisson")  pq = 1.0/mu;

    NumericVector sn;
    if(q_Nested > 0){
      IntegerVector idx2 = wrap(seq(q_nn, q_nn + q_Nested - 1));
      NumericVector sn0 = par[idx2];
      rowvec sn1 = real(as<rowvec>(sn0));
      sn = as<NumericVector>(wrap(sn1));
    }

    // Sparse A → dense (identical path to pglmm_iV_logdetV_cpp)
    arma::sp_mat A_sp;
    if(family=="binomial") A_sp = sp_mat(diagmat(pq));
    if(family=="poisson")  A_sp = sp_mat(diagmat(1.0/mu));
    for(int j = 0; j < q_Nested; ++j){
      double snj2 = pow((double)sn[j], 2);
      sp_mat nj = nested[j];
      A_sp = A_sp + snj2 * nj;
    }
    arma::mat A1(A_sp);

    arma::mat chol_A;
    if(!arma::chol(chol_A, A1)){
      return 1e15;   // not PD; optimizer will step away
    }
    arma::mat Rlo = chol_A.t();   // lower L s.t. L*L' = A
    double logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));

    // iA*H via two triangular solves (O(n^2), no explicit inversion)
    arma::vec iA_H = arma::solve(arma::trimatl(Rlo), H);
    iA_H = arma::solve(arma::trimatu(chol_A), iA_H);

    double logdetV;
    arma::vec iVH;

    if(q_nn > 0){
      arma::sp_mat Ut_sp = build_Ut_helper(par, St, Zt, q_nn);
      arma::mat Ut_d = arma::mat(Ut_sp);
      arma::mat U_d  = Ut_d.t();
      int Q = (int)Ut_d.n_rows;

      // iA*U
      arma::mat iA_U = arma::solve(arma::trimatl(Rlo), U_d);
      iA_U = arma::solve(arma::trimatu(chol_A), iA_U);

      // M = (I + Ut*iA*U)^{-1}
      arma::mat IpUtAU = arma::eye(Q, Q) + Ut_d * iA_U;
      arma::mat M      = arma::inv_sympd(IpUtAU);

      // iVH = iA_H - iA_U*(M*(Ut*iA_H))   [Woodbury, no n×n iV]
      iVH = iA_H - iA_U * (M * (Ut_d * iA_H));

      // logdetV = logdetA + log|I + Ut*iA*U|  (matrix-det lemma)
      double sg; arma::log_det(logdetV, sg, IpUtAU);
      logdetV += logdetA;

      if(REML){
        arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
        iA_X = arma::solve(arma::trimatu(chol_A), iA_X);
        arma::mat iVX = iA_X - iA_U * (M * (Ut_d * iA_X));
        double ld, sg2; arma::log_det(ld, sg2, X.t() * iVX);
        LL = 0.5*(logdetV + arma::dot(H, iVH) + ld);
      } else {
        LL = 0.5*(logdetV + arma::dot(H, iVH));
      }
    } else {
      // q_nn == 0: iV = iA
      iVH     = iA_H;
      logdetV = logdetA;
      if(REML){
        arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
        iA_X = arma::solve(arma::trimatu(chol_A), iA_X);
        double ld, sg; arma::log_det(ld, sg, X.t() * iA_X);
        LL = 0.5*(logdetV + arma::dot(H, iVH) + ld);
      } else {
        LL = 0.5*(logdetV + arma::dot(H, iVH));
      }
    }
  }

  if(verbose) Rcout << LL << " " << par << std::endl;
  return LL;
}

// ============================================================
// OPT-K: PQL analytic gradient for L-BFGS
//
// The PQL inner step minimises:
//   LL_PQL = 0.5*(logdetV + H'iVH + [REML: logdet(X'iVX)])
// where V = diag(iW) + random-effects-structure(par),  iW = 1/W = 1/(size*mu*(1-mu)).
// H is the pre-computed pseudo-residual, fixed for the duration of this optimizer call.
//
// Gradient differs from gauss_ll_and_grad in three ways:
//   (1) Residual diagonal is diag(iW), not s2resid*I — heteroscedastic.
//   (2) Quadratic coefficient is 1.0, not factor/HiVH — no LL concentration over s2resid.
//   (3) H is supplied; B is not re-estimated.
//
// Key identity (holds for any A s.t. V = A + Ut'*Ut):
//   trace_diag_k = u_k'*iV*u_k = (Ut*A^{-1}*U).diag()[k]
//                               - diag((Ut*A^{-1}*U)' * M * (Ut*A^{-1}*U))[k]
//   where M = inv(I + Ut*A^{-1}*U).
// For Gaussian A=I this reduces to 1−M_kk; verified analytically.
// ============================================================

struct PqlCtx {
    const arma::vec*    H;      // pre-computed pseudo-residual (fixed)
    const arma::mat*    X;      // fixed-effects design matrix (for REML)
    const arma::sp_mat* Zt;
    const arma::sp_mat* St;
    const Rcpp::List*   nested;
    const arma::vec*    W;      // working weights: totalSize*mu*(1-mu) or mu
    bool REML;
    bool verbose;
};

static double pql_ll_and_grad(unsigned n_par, const double* x, double* grad,
                               const PqlCtx& d) {
    const arma::vec&    H      = *d.H;
    const arma::mat&    X      = *d.X;
    const arma::sp_mat& Zt     = *d.Zt;
    const arma::sp_mat& St     = *d.St;
    const Rcpp::List&   nested = *d.nested;
    const arma::vec&    W      = *d.W;
    bool REML = d.REML;

    int n        = (int)H.n_elem;
    int p        = (int)X.n_cols;
    int q_nn     = (int)St.n_rows;
    int q_Nested = (int)nested.size();
    int q_total  = q_nn + q_Nested;

    for (int k = 0; k < q_total; ++k) grad[k] = 0.0;

    // extract parameters (abs: optimizer may pass negative values)
    arma::vec iC1;
    arma::mat Ut_d, U_d, St_d;
    int Q = 0;
    if (q_nn > 0) {
        arma::rowvec sr(q_nn);
        for (int r = 0; r < q_nn; ++r) sr[r] = std::abs(x[r]);
        arma::mat iC0 = sr * arma::mat(St);
        iC1 = arma::vectorise(iC0, 0);
        Q   = (int)iC1.n_elem;
        arma::sp_mat Ut_sp = arma::sp_mat(arma::diagmat(iC1)) * Zt;
        Ut_d = arma::mat(Ut_sp);
        U_d  = Ut_d.t();
        St_d = arma::mat(St);
    }
    arma::vec sn_vec(q_Nested);
    for (int j = 0; j < q_Nested; ++j)
        sn_vec[j] = std::abs(x[q_nn + j]);

    // block detection for nested (used for Path B / B+A)
    int n_sp = 0, n_site = 0;
    bool use_blk = false;
    if (q_Nested > 0) {
        arma::sp_mat nj0_sp = nested[0];
        n_sp   = find_block_size(nj0_sp, n);
        n_site = (n_sp > 0) ? n / n_sp : 0;
        use_blk = (n_sp > 0);
    }

    double logdetV = 0.0, HiVH = 0.0;
    arma::vec v;
    arma::mat denom, iVX_d;
    bool have_reml_denom = false;

    // Woodbury state (non-nested gradient)
    arma::mat M_mat, ZtiAU_g;
    bool have_M = false;

    // block-Chol state (nested gradient)
    std::vector<arma::mat> chols_blk;

    // dense Chol / LU state
    arma::mat chol_A_dense;
    bool have_chol_dense = false;
    arma::mat iA_U_g;
    arma::mat iV_dense;
    bool have_iV_dense = false;

    // ===== PATH A: non-nested only, Woodbury with D=diag(iW) =====
    if (q_Nested == 0) {
        // WU_d = U_d with each column scaled by W  (analogous to iA*U when A=diag(iW))
        arma::mat WU_d = U_d;
        WU_d.each_col() %= W;

        arma::mat IpUtWU = arma::symmatu(arma::eye(Q, Q) + Ut_d * WU_d);
        arma::mat M = arma::inv_sympd(IpUtWU);

        double lsign;
        arma::log_det(logdetV, lsign, IpUtWU);
        if (!std::isfinite(logdetV)) {
            arma::mat lg = arma::chol(IpUtWU);
            logdetV = 2.0 * arma::accu(arma::log(lg.diag()));
        }
        logdetV -= arma::accu(arma::log(W));   // logdet(diag(W)^{-1}) = -sum(log W)

        arma::vec WH = W % H;
        v    = WH - WU_d * (M * (Ut_d * WH));
        HiVH = arma::dot(H, v);

        if (REML) {
            arma::mat WX = X;
            WX.each_col() %= W;
            iVX_d = WX - WU_d * (M * (Ut_d * WX));
            denom = X.t() * iVX_d;
            have_reml_denom = true;
        }

        M_mat   = M;
        ZtiAU_g = IpUtWU - arma::eye(Q, Q);   // = Ut_d * WU_d = Ut * A^{-1} * U
        have_M  = true;

    // ===== PATH B: nested-only, block Cholesky with pq=iW =====
    } else if (q_nn == 0 && use_blk) {
        std::vector<arma::mat> nj_fixed(n_site,
                                        arma::mat(n_sp, n_sp, arma::fill::zeros));
        for (int j = 0; j < q_Nested; ++j) {
            double snj2 = sn_vec[j] * sn_vec[j];
            arma::sp_mat njsp = nested[j];
            for (int k = 0; k < n_site; ++k) {
                int s = k * n_sp;
                nj_fixed[k] += snj2 *
                    arma::mat(njsp(arma::span(s, s+n_sp-1), arma::span(s, s+n_sp-1)));
            }
        }
        chols_blk.resize(n_site);
        arma::vec pq = 1.0 / W;    // heteroscedastic residual diagonal
        double logdetA = block_chol(chols_blk, pq, nj_fixed, n_sp, n_site);
        logdetV = logdetA;

        v.set_size(n);
        for (int k = 0; k < n_site; ++k) {
            int s = k * n_sp;
            arma::mat Rlo_k = chols_blk[k].t();
            arma::vec tmp = arma::solve(arma::trimatl(Rlo_k), H.subvec(s, s+n_sp-1));
            v.subvec(s, s+n_sp-1) = arma::solve(arma::trimatu(chols_blk[k]), tmp);
        }
        HiVH = arma::dot(H, v);

        if (REML) {
            iVX_d.set_size(n, p);
            for (int k = 0; k < n_site; ++k) {
                int s = k * n_sp;
                arma::mat Rlo_k = chols_blk[k].t();
                arma::mat tmp = arma::solve(arma::trimatl(Rlo_k), X.rows(s, s+n_sp-1));
                iVX_d.rows(s, s+n_sp-1) = arma::solve(arma::trimatu(chols_blk[k]), tmp);
            }
            denom = X.t() * iVX_d;
            have_reml_denom = true;
        }

    // ===== PATH B+A: mixed nested+non-nested, block Cholesky + Woodbury =====
    // A = diag(1/W) + sum_j sn_j² * nj is block-diagonal in site-major order.
    // Block-Chol A is O(n_site*n_sp³); Woodbury handles the non-nested columns.
    } else if (use_blk) {
        std::vector<arma::mat> nj_fixed(n_site,
                                        arma::mat(n_sp, n_sp, arma::fill::zeros));
        for (int j = 0; j < q_Nested; ++j) {
            double snj2 = sn_vec[j] * sn_vec[j];
            arma::sp_mat njsp = nested[j];
            for (int k = 0; k < n_site; ++k) {
                int s = k * n_sp;
                nj_fixed[k] += snj2 *
                    arma::mat(njsp(arma::span(s, s+n_sp-1), arma::span(s, s+n_sp-1)));
            }
        }
        chols_blk.resize(n_site);
        arma::vec pq = 1.0 / W;
        double logdetA = block_chol(chols_blk, pq, nj_fixed, n_sp, n_site);

        // iA_H and iA_U via block triangular solves
        arma::vec iA_H(n);
        arma::mat iA_U(n, Q);
        for (int k = 0; k < n_site; ++k) {
            int s = k * n_sp;
            arma::mat Rlo_k = chols_blk[k].t();
            arma::vec hk = arma::solve(arma::trimatl(Rlo_k), H.subvec(s, s+n_sp-1));
            iA_H.subvec(s, s+n_sp-1) = arma::solve(arma::trimatu(chols_blk[k]), hk);
            arma::mat uk = arma::solve(arma::trimatl(Rlo_k), U_d.rows(s, s+n_sp-1));
            iA_U.rows(s, s+n_sp-1) = arma::solve(arma::trimatu(chols_blk[k]), uk);
        }
        iA_U_g = iA_U;

        arma::mat ZtiAU = Ut_d * iA_U;
        ZtiAU_g = ZtiAU;
        arma::mat IpUtAU = arma::symmatu(arma::eye(Q, Q) + ZtiAU);
        arma::mat M = arma::inv_sympd(IpUtAU);
        M_mat = M; have_M = true;

        double sv; arma::log_det(logdetV, sv, IpUtAU);
        logdetV += logdetA;

        v = iA_H - iA_U * (M * (Ut_d * iA_H));
        HiVH = arma::dot(H, v);

        if (REML) {
            arma::mat iA_X(n, p);
            for (int k = 0; k < n_site; ++k) {
                int s = k * n_sp;
                arma::mat Rlo_k = chols_blk[k].t();
                arma::mat xk = arma::solve(arma::trimatl(Rlo_k), X.rows(s, s+n_sp-1));
                iA_X.rows(s, s+n_sp-1) = arma::solve(arma::trimatu(chols_blk[k]), xk);
            }
            iVX_d = iA_X - iA_U * (M * (Ut_d * iA_X));
            denom = X.t() * iVX_d;
            have_reml_denom = true;
        }

    // ===== PATH C: full dense (non-block nested or fallback) =====
    } else {
        arma::vec pq = 1.0 / W;
        arma::sp_mat A_sp = arma::sp_mat(arma::diagmat(pq));
        for (int j = 0; j < q_Nested; ++j)
            A_sp = A_sp + sn_vec[j]*sn_vec[j] * arma::sp_mat(nested[j]);
        arma::mat A1(A_sp);

        arma::mat chol_A;
        if (!arma::chol(chol_A, A1)) {
            // LU fallback
            arma::mat iA1 = arma::inv(A1);
            double logdetA, sgn; arma::log_det(logdetA, sgn, A1);
            if (q_nn > 0) {
                arma::mat iA_U = iA1 * U_d;
                arma::mat ZtiAU = Ut_d * iA_U;
                arma::mat Ishort = arma::symmatu(arma::eye(Q, Q) + ZtiAU);
                arma::mat M = arma::inv_sympd(Ishort);
                M_mat = M; ZtiAU_g = ZtiAU; iA_U_g = iA_U; have_M = true;
                iV_dense = iA1 - iA_U * M * iA_U.t();
                double sv; arma::log_det(logdetV, sv, Ishort); logdetV += logdetA;
            } else {
                iV_dense = iA1; logdetV = logdetA;
            }
            have_iV_dense = true;
            v    = iV_dense * H;
            HiVH = arma::dot(H, v);
            if (REML) { iVX_d = iV_dense * X; denom = X.t() * iVX_d; have_reml_denom = true; }
        } else {
            chol_A_dense = chol_A; have_chol_dense = true;
            arma::mat Rlo = chol_A.t();
            double logdetA = 2.0 * arma::accu(arma::log(chol_A.diag()));
            if (q_nn > 0) {
                arma::mat iA_U = arma::solve(arma::trimatl(Rlo), U_d);
                iA_U           = arma::solve(arma::trimatu(chol_A), iA_U);
                iA_U_g         = iA_U;
                arma::mat ZtiAU = Ut_d * iA_U;
                ZtiAU_g = ZtiAU;
                arma::mat Ishort = arma::symmatu(arma::eye(Q, Q) + ZtiAU);
                arma::mat M = arma::inv_sympd(Ishort);
                M_mat = M; have_M = true;
                arma::mat iA_H = arma::solve(arma::trimatl(Rlo), H);
                iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
                arma::vec Ut_iA_H = Ut_d * iA_H;
                v    = iA_H - iA_U * (M * Ut_iA_H);
                HiVH = arma::dot(H, iA_H) - arma::as_scalar(Ut_iA_H.t() * (M * Ut_iA_H));
                double sv; arma::log_det(logdetV, sv, Ishort); logdetV += logdetA;
                if (REML) {
                    arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
                    iA_X           = arma::solve(arma::trimatu(chol_A), iA_X);
                    iVX_d = iA_X - iA_U * (M * (Ut_d * iA_X));
                    denom = X.t() * iVX_d; have_reml_denom = true;
                }
            } else {
                arma::mat iA_H = arma::solve(arma::trimatl(Rlo), H);
                iA_H           = arma::solve(arma::trimatu(chol_A), iA_H);
                v = iA_H; HiVH = arma::dot(H, v); logdetV = logdetA;
                if (REML) {
                    arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
                    iA_X           = arma::solve(arma::trimatu(chol_A), iA_X);
                    iVX_d = iA_X; denom = X.t() * iVX_d; have_reml_denom = true;
                }
            }
        }
    }

    // ===== LL =====
    double LL = 0.5 * (logdetV + HiVH);
    if (REML && have_reml_denom) {
        double ld, sg; arma::log_det(ld, sg, denom);
        LL += 0.5 * ld;
    }
    if (!std::isfinite(LL)) return 1e30;

    // ===== REML F matrices (for gradient) =====
    arma::mat F_mat;
    arma::vec rowNorm2_UtF(Q > 0 ? Q : 1, arma::fill::zeros);
    bool F_built = false;
    if (REML && q_nn > 0 && have_reml_denom) {
        arma::mat denom_inv = arma::solve(denom, arma::eye(p, p));
        F_mat = iVX_d * denom_inv;
        arma::mat UtF = Ut_d * F_mat;
        rowNorm2_UtF = arma::sum((UtF * denom) % UtF, 1);
        F_built = true;
    }
    arma::mat F_mat_nested;
    bool F_nested_built = false;
    if (REML && q_Nested > 0 && have_reml_denom) {
        if (F_built) {
            F_mat_nested = F_mat;
        } else {
            arma::mat denom_inv = arma::solve(denom, arma::eye(p, p));
            F_mat_nested = iVX_d * denom_inv;
        }
        F_nested_built = true;
    }

    // ===== Gradient: non-nested sr_r =====
    if (q_nn > 0 && have_M) {
        // trace_diag_pql[k] = (Ut*A^{-1}*U).diag()[k] - diag((Ut*A^{-1}*U)' * M * (Ut*A^{-1}*U))[k]
        // ZtiAU_g = Ut * A^{-1} * U  (Q×Q, computed in each path above)
        arma::mat MC = M_mat * ZtiAU_g;
        arma::vec diag_CtMC = arma::sum(ZtiAU_g % MC, 0).t();
        arma::vec trace_diag = ZtiAU_g.diag() - diag_CtMC;

        arma::vec b = Ut_d * v;   // Q-vector: Ut * iV * H

        arma::vec tc_raw = trace_diag - (b % b) -
                           (F_built ? rowNorm2_UtF : arma::zeros<arma::vec>(Q));

        arma::vec tc_scaled(Q);
        for (int k = 0; k < Q; ++k)
            tc_scaled[k] = (std::abs(iC1[k]) > 1e-15) ? tc_raw[k] / iC1[k] : 0.0;

        for (int r = 0; r < q_nn; ++r)
            grad[r] = arma::dot(St_d.row(r), tc_scaled);
    }

    // ===== Gradient: nested sn_j =====
    if (q_Nested > 0) {
        for (int j = 0; j < q_Nested; ++j) {
            arma::sp_mat nj_sp = nested[j];
            double tr_iV_nj = 0.0;

            if (!chols_blk.empty()) {
                // PATH B / B+A: tr(iA * nj) via block inverses
                for (int k = 0; k < n_site; ++k) {
                    int s = k * n_sp;
                    arma::mat nj_blk = arma::mat(
                        nj_sp(arma::span(s, s+n_sp-1), arma::span(s, s+n_sp-1)));
                    tr_iV_nj += arma::accu(phyr_chol2inv(chols_blk[k]) % nj_blk);
                }
                // PATH B+A only: Woodbury correction tr(iA_U' * nj * iA_U * M)
                if (have_M) {
                    arma::mat NjiAU = nj_sp * iA_U_g;
                    tr_iV_nj -= arma::accu((iA_U_g.t() * NjiAU) % M_mat);
                }
            } else if (have_chol_dense) {
                arma::mat Rlo = chol_A_dense.t();
                arma::mat Nj_d = arma::mat(nj_sp);
                arma::mat W_s = arma::solve(arma::trimatl(Rlo), Nj_d);
                W_s = arma::solve(arma::trimatu(chol_A_dense), W_s);
                double tr_iA_nj = arma::trace(W_s);
                if (q_nn > 0 && have_M) {
                    arma::mat NjiAU = nj_sp * iA_U_g;
                    double tr_corr = arma::accu((iA_U_g.t() * NjiAU) % M_mat);
                    tr_iV_nj = tr_iA_nj - tr_corr;
                } else {
                    tr_iV_nj = tr_iA_nj;
                }
            } else {
                tr_iV_nj = arma::accu(iV_dense % arma::mat(nj_sp));
            }

            double vtnjv = arma::dot(v, nj_sp * v);
            double reml_nj = 0.0;
            if (REML && F_nested_built)
                reml_nj = arma::accu(iVX_d % (nj_sp * F_mat_nested));

            // PQL coefficient is 1.0 (not factor/HiVH as in Gaussian concentrated LL)
            grad[q_nn + j] = sn_vec[j] * (tr_iV_nj - vtnjv - reml_nj);
        }
    }

    for (int k = 0; k < q_total; ++k)
        if (!std::isfinite(grad[k])) grad[k] = 0.0;

    if (d.verbose) Rcpp::Rcout << "pql_grad LL=" << LL << std::endl;
    return LL;
}

static double pql_ll_nlopt_cb(unsigned n_par, const double* x,
                               double* grad, void* f_data) {
    const PqlCtx* d = static_cast<const PqlCtx*>(f_data);
    if (grad) {
        try { return pql_ll_and_grad(n_par, x, grad, *d); }
        catch (...) { for (unsigned k = 0; k < n_par; ++k) grad[k] = 0.0; return 1e30; }
    }
    try {
        std::vector<double> tmp(n_par, 0.0);
        return pql_ll_and_grad(n_par, x, tmp.data(), *d);
    } catch (...) { return 1e30; }
}

// ============================================================
// pglmm_internal_cpp: OPT5b + OPT6b inner loops
// ============================================================
// [[Rcpp::export]]
List pglmm_internal_cpp(const arma::mat& X, const arma::vec& Y,
                        const arma::sp_mat& Zt, const arma::sp_mat& St,
                        const List& nested, const bool REML, const bool verbose,
                        const int n, const int p, const int q, const int maxit,
                        const double reltol, const double tol_pql,
                        const double maxit_pql, const std::string optimizer,
                        arma::mat B_init, arma::vec ss,
                        const std::string family, arma::vec totalSize){
  Rcpp::checkUserInterrupt();
  mat B = B_init;
  mat b(n, 1, fill::zeros);
  mat beta = join_vert(B, b);
  vec mu;
  if(family=="binomial") mu = arma::exp(X*B)/(1+arma::exp(X*B));
  if(family=="poisson")  mu = arma::exp(X*B);
  mat XX = join_horiz(X, arma::eye(n,n));

  vec est_ss = ss, est_B = B;
  vec oldest_ss(ss.n_elem); oldest_ss.fill(1e6);
  mat oldest_B(B.n_rows, B.n_cols); oldest_B.fill(1e6);

  int q_nn     = St.n_rows;
  int q_Nested = (int)nested.size();

  // Detect block structure once (geometry doesn't change).
  // use_blk is only true when data are sorted site×species AND the nested
  // matrix is genuinely block-diagonal; verify_block_structure guards against
  // silently wrong results when data ordering doesn't match expectations.
  int n_sp = 0, n_site = 0;
  bool use_blk = false;
  if(q_Nested > 0){
    sp_mat nj0_sp = nested[0];   // implicit conversion
    n_sp   = find_block_size(nj0_sp, n);
    n_site = (n_sp > 0) ? n / n_sp : 0;
    use_blk = (n_sp > 0);
  }

  unsigned int iteration = 0;
  double tol2 = tol_pql * tol_pql;
  double LL   = 0.0;
  NumericVector ss0 = wrap(ss);
  vec Z, H, niter;
  int convcode = 0;
  mat iV;

  Rcpp::Environment stats_env("package:stats");
  Rcpp::Function optim_fn = stats_env["optim"];
  Rcpp::Environment nloptr_pkg = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptr_pkg["nloptr"];
  Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
  Rcpp::Function pglmm_LL_cpp2 = phyr_pkg["pglmm_LL_cpp"];

  // ==== outer PQL loop ====
  while((as_scalar(trans(est_ss-oldest_ss)*(est_ss-oldest_ss)) > tol2 ||
         as_scalar(trans(est_B -oldest_B )*(est_B -oldest_B )) > tol2) &&
        iteration <= maxit_pql){
    oldest_ss = est_ss; oldest_B = est_B;
    Rcpp::checkUserInterrupt();

    unsigned int iteration_m = 0;
    vec est_B_m = B;
    mat oldest_B_m(B.n_rows, B.n_cols); oldest_B_m.fill(1e6);

    // ==== inner mean loop — choose path based on structure ====

    if(q_Nested == 0 && q_nn > 0){
      // OPT5b: diagonal iA, never form n×n iV
      arma::sp_mat Ut_sp = build_Ut_helper(ss0, St, Zt, q_nn);
      arma::mat Ut_d = arma::mat(Ut_sp);
      arma::mat U_d  = Ut_d.t();
      int Q = (int)Ut_d.n_rows;

      while(as_scalar(trans(est_B_m-oldest_B_m)*(est_B_m-oldest_B_m)) > tol2 &&
            iteration_m <= maxit_pql){
        oldest_B_m = est_B_m;

        arma::vec W;
        if(family=="binomial") W = totalSize % mu % (1.0-mu);
        if(family=="poisson")  W = mu;

        arma::mat WU = U_d; WU.each_col() %= W;
        arma::mat M  = arma::inv_sympd(arma::eye(Q,Q) + Ut_d * WU);

        if(family=="binomial") Z = X*B + b + (Y/totalSize-mu)/(mu%(1.0-mu));
        if(family=="poisson")  Z = X*B + b + (Y-mu)/mu;

        arma::mat WX = X; WX.each_col() %= W;
        arma::mat iVX   = WX  - WU * (M * (Ut_d * WX));
        arma::mat denom = X.t() * iVX;

        arma::vec WZ  = W % Z;
        arma::vec iVZ = WZ - WU * (M * (Ut_d * WZ));
        arma::mat num = X.t() * iVZ;

        B = solve(denom, num);

        arma::vec r   = Z - X * B;
        arma::vec Wr  = W % r;
        arma::vec iVr = Wr - WU * (M * (Ut_d * Wr));
        b = U_d * (Ut_d * iVr);

        beta = join_vert(B, b);
        if(family=="binomial") mu = arma::exp(XX*beta)/(1+arma::exp(XX*beta));
        if(family=="poisson")  mu = arma::exp(XX*beta);
        est_B_m = B;
        if(verbose) Rcout << "mean part: " << iteration_m << " " << trans(B);
        ++iteration_m;
        if(B.has_nan()) Rcpp::stop("Estimation of B failed.");
      }

    } else if(q_Nested > 0 && use_blk && q_nn == 0){
      // OPT6b: nested-only, block-diagonal A → block Cholesky (O(n_site*n_sp^3) vs O(n^3)).
      // Restricted to q_nn==0: when q_nn>0 the Woodbury step amplifies tiny iA
      // differences between block-Chol and full-Chol, causing ss divergence vs fallback.

      // sn from current ss0 (fixed this outer iteration)
      arma::vec sn_vec(q_Nested);
      for(int j = 0; j < q_Nested; ++j)
        sn_vec[j] = std::abs((double)ss0[q_nn + j]);

      // Precompute nested fixed part per block (changes only when ss0 changes)
      std::vector<arma::mat> nj_fixed(n_site,
                                      arma::mat(n_sp, n_sp, arma::fill::zeros));
      for(int j = 0; j < q_Nested; ++j){
        double snj2 = sn_vec[j] * sn_vec[j];
        sp_mat njsp = nested[j];   // implicit conversion
        for(int k = 0; k < n_site; ++k){
          int s = k * n_sp;
          nj_fixed[k] += snj2 *
            arma::mat(njsp(arma::span(s, s+n_sp-1), arma::span(s, s+n_sp-1)));
        }
      }

      // Non-nested Ut/U (fixed this outer iteration)
      arma::mat Ut_d, U_d;
      int Q = 0;
      if(q_nn > 0){
        arma::sp_mat Ut_sp = build_Ut_helper(ss0, St, Zt, q_nn);
        Ut_d = arma::mat(Ut_sp); U_d = Ut_d.t();
        Q = (int)Ut_d.n_rows;
      }

      while(as_scalar(trans(est_B_m-oldest_B_m)*(est_B_m-oldest_B_m)) > tol2 &&
            iteration_m <= maxit_pql){
        oldest_B_m = est_B_m;

        // pq changes with mu each inner step
        arma::vec pq;
        if(family=="binomial") pq = 1.0/(totalSize % mu % (1.0-mu));
        if(family=="poisson")  pq = 1.0/mu;

        // Block Cholesky: O(n_site × n_sp³) — speedup vs O(n³) full dense Chol
        std::vector<arma::mat> chols(n_site);
        block_chol(chols, pq, nj_fixed, n_sp, n_site);

        // Assemble full dense block-diagonal iA from block Cholesky inverses
        arma::mat iA_dense(n, n, arma::fill::zeros);
        for(int k = 0; k < n_site; ++k){
          int s = k * n_sp;
          arma::mat Ri = arma::inv(arma::trimatu(chols[k]));
          iA_dense(arma::span(s, s+n_sp-1), arma::span(s, s+n_sp-1)) = Ri * Ri.t();
        }

        // Woodbury → dense iV, replicating pglmm_iV_logdetV_cpp sparse path exactly
        arma::mat iV_dense;
        if(Q > 0){
          // Use sparse iA (matches sp_mat(Ri*Ri.t()) in pglmm_iV_logdetV_cpp)
          arma::sp_mat iA_sp = arma::sp_mat(iA_dense);
          arma::sp_mat Ut_sp_loc = build_Ut_helper(ss0, St, Zt, q_nn);
          arma::sp_mat U_sp_loc  = arma::trans(Ut_sp_loc);
          arma::sp_mat Ish_sp((arma::uword)Q, (arma::uword)Q); Ish_sp.eye();
          arma::mat Ish_UtAU = arma::mat(Ish_sp + Ut_sp_loc * iA_sp * U_sp_loc);
          arma::sp_mat M_sp  = arma::sp_mat(arma::inv_sympd(Ish_UtAU));
          arma::sp_mat iV0   = iA_sp - iA_sp * U_sp_loc * M_sp * Ut_sp_loc * iA_sp;
          iV_dense = arma::mat(iV0);
        } else {
          iV_dense = iA_dense;
        }

        if(family=="binomial") Z = X*B + b + (Y/totalSize-mu)/(mu%(1.0-mu));
        if(family=="poisson")  Z = X*B + b + (Y-mu)/mu;

        arma::mat denom = X.t() * (iV_dense * X);
        arma::mat num   = X.t() * (iV_dense * Z);
        B = solve(denom, num);

        arma::vec r    = Z - X * B;
        arma::vec iV_r = iV_dense * r;

        // b = (V - diag(pq)) * iV_r = sum_j snj^2 * nj * iV_r + U*Ut*iV_r
        arma::vec b_vec(n, arma::fill::zeros);
        for(int j = 0; j < q_Nested; ++j){
          double snj2 = sn_vec[j] * sn_vec[j];
          sp_mat njsp = nested[j];   // implicit conversion
          b_vec += snj2 * (njsp * iV_r);
        }
        if(Q > 0) b_vec += U_d * (Ut_d * iV_r);
        b = b_vec;

        beta = join_vert(B, b);
        if(family=="binomial") mu = arma::exp(XX*beta)/(1+arma::exp(XX*beta));
        if(family=="poisson")  mu = arma::exp(XX*beta);
        est_B_m = B;
        if(verbose) Rcout << "mean part: " << iteration_m << " " << trans(B);
        ++iteration_m;
        if(B.has_nan()) Rcpp::stop("Estimation of B failed.");
      }

    } else {
      // Fallback: triangular-solve path (same A construction as pglmm_iV_logdetV_cpp,
      // but triangular solves instead of inv()+sp_mat — no n×n iV materialized).

      // Extract sn (nested variances) from ss0 — fixed this outer iteration
      arma::vec sn_fb(q_Nested);
      {
        IntegerVector idx2 = wrap(seq(q_nn, q_nn + q_Nested - 1));
        NumericVector sn0 = ss0[idx2];
        rowvec sn1 = real(as<rowvec>(sn0));
        sn_fb = as<arma::vec>(wrap(sn1));
      }
      // Non-nested Ut/U from ss0
      arma::mat Ut_fb, U_fb;
      int Q_fb = 0;
      if(q_nn > 0){
        arma::sp_mat Ut_sp = build_Ut_helper(ss0, St, Zt, q_nn);
        Ut_fb = arma::mat(Ut_sp); U_fb = Ut_fb.t();
        Q_fb = (int)Ut_fb.n_rows;
      }

      while(as_scalar(trans(est_B_m-oldest_B_m)*(est_B_m-oldest_B_m)) > tol2 &&
            iteration_m <= maxit_pql){
        oldest_B_m = est_B_m;

        // Build A sparse → dense (same path as pglmm_iV_logdetV_cpp)
        arma::sp_mat A_sp;
        if(family=="binomial") A_sp = sp_mat(diagmat(1.0/(totalSize%mu%(1.0-mu))));
        if(family=="poisson")  A_sp = sp_mat(diagmat(1.0/mu));
        for(int j = 0; j < q_Nested; ++j){
          double snj2 = sn_fb[j] * sn_fb[j];
          sp_mat nj = nested[j];
          A_sp = A_sp + snj2 * nj;
        }
        arma::mat A1(A_sp);
        arma::mat chol_A;
        if (!arma::chol(chol_A, A1)) break;  // A not PD; use last valid estimates
        arma::mat Rlo = chol_A.t();

        if(family=="binomial") Z = X*B + b + (Y/totalSize-mu)/(mu%(1.0-mu));
        if(family=="poisson")  Z = X*B + b + (Y-mu)/mu;

        // iA*X and iA*Z via triangular solves
        arma::mat iA_X = arma::solve(arma::trimatl(Rlo), X);
        iA_X = arma::solve(arma::trimatu(chol_A), iA_X);
        arma::vec iA_Z = arma::solve(arma::trimatl(Rlo), Z);
        iA_Z = arma::solve(arma::trimatu(chol_A), iA_Z);

        arma::mat denom, iA_U_fb;
        arma::mat M_fb;
        if(Q_fb > 0){
          iA_U_fb = arma::solve(arma::trimatl(Rlo), U_fb);
          iA_U_fb = arma::solve(arma::trimatu(chol_A), iA_U_fb);
          arma::mat IpUtAU = arma::eye(Q_fb,Q_fb) + Ut_fb * iA_U_fb;
          M_fb  = arma::inv_sympd(IpUtAU);
          arma::mat iVX = iA_X - iA_U_fb*(M_fb*(Ut_fb*iA_X));
          arma::vec iVZ = iA_Z - iA_U_fb*(M_fb*(Ut_fb*iA_Z));
          denom = X.t() * iVX;
          arma::mat num = X.t() * iVZ;
          B = arma::solve(denom, num);
        } else {
          denom = X.t() * iA_X;
          arma::mat num = X.t() * iA_Z;
          B = arma::solve(denom, num);
        }

        arma::vec r = Z - X * B;
        arma::vec iA_r = arma::solve(arma::trimatl(Rlo), r);
        iA_r = arma::solve(arma::trimatu(chol_A), iA_r);
        arma::vec iV_r = (Q_fb > 0) ? iA_r - iA_U_fb*(M_fb*(Ut_fb*iA_r)) : iA_r;

        // b = (V - iW)*iV_r = sum_j snj^2 * nj * iV_r + U*Ut*iV_r
        arma::vec b_vec(n, arma::fill::zeros);
        for(int j = 0; j < q_Nested; ++j){
          double snj2 = sn_fb[j] * sn_fb[j];
          sp_mat nj = nested[j];
          b_vec += snj2 * (nj * iV_r);
        }
        if(Q_fb > 0) b_vec += U_fb * (Ut_fb * iV_r);
        b = b_vec;

        beta = join_vert(B, b);
        if(family=="binomial") mu = arma::exp(XX*beta)/(1+arma::exp(XX*beta));
        if(family=="poisson")  mu = arma::exp(XX*beta);
        est_B_m = B;
        if(verbose) Rcout << "mean part: " << iteration_m << " " << trans(B);
        ++iteration_m;
        if(B.has_nan()) Rcpp::stop("Estimation of B failed.");
      }
    }

    // ==== variance component (outer optimiser) ====
    if(family=="binomial") Z = X*B + b + (Y/totalSize-mu)/(mu%(1-mu));
    if(family=="poisson")  Z = X*B + b + (Y-mu)/mu;
    H = Z - X * B;

    Rcpp::List opt;
    if(optimizer == "Nelder-Mead"){
      if(q > 1){
        opt = optim_fn(_["par"]=ss0, _["fn"]=pglmm_LL_cpp2,
                       _["H"]=H, _["X"]=X, _["Zt"]=Zt, _["St"]=St,
                       _["mu"]=mu, _["nested"]=nested,
                       _["REML"]=REML, _["verbose"]=verbose,
                       _["family"]=family, _["totalSize"]=totalSize,
                       _["method"]="Nelder-Mead",
                       _["control"]=List::create(_["maxit"]=maxit,_["reltol"]=reltol));
      } else {
        Rcpp::stop("With 1 random term and cpp=TRUE use a different optimizer.");
      }
    } else if (optimizer == "lbfgs") {
      // OPT-K: NLopt C API with analytic PQL gradient
      arma::vec W_vec;
      if(family=="binomial") W_vec = totalSize % mu % (1.0 - mu);
      if(family=="poisson")  W_vec = mu;
      PqlCtx ctx = {&H, &X, &Zt, &St, &nested, &W_vec, REML, verbose};
      nlopt_opt opt_c = nlopt_create(NLOPT_LD_LBFGS, (unsigned)q);
      nlopt_set_min_objective(opt_c, pql_ll_nlopt_cb, &ctx);
      // Variance parameters are non-negative; bounding [0, ∞) prevents lbfgs
      // from crossing zero where abs() creates a gradient sign flip that
      // causes PQL outer-loop oscillation.
      std::vector<double> lb_vec(q, 0.0);
      nlopt_set_lower_bounds(opt_c, lb_vec.data());
      nlopt_set_ftol_rel(opt_c, reltol);
      nlopt_set_ftol_abs(opt_c, reltol);
      nlopt_set_xtol_rel(opt_c, 0.0001);
      nlopt_set_maxeval(opt_c, maxit);
      // Starting point: ensure strictly positive (gradient undefined at 0)
      std::vector<double> x_vec(q);
      for (int k = 0; k < q; ++k) x_vec[k] = std::max((double)ss0[k], 1e-6);
      double opt_val = 1e30;
      nlopt_result res = nlopt_optimize(opt_c, x_vec.data(), &opt_val);
      int cc = (int)res;
      nlopt_destroy(opt_c);
      arma::vec par0_lbfgs = arma::abs(arma::vec(x_vec.data(), (arma::uword)q));
      opt = List::create(_["par"]=Rcpp::wrap(par0_lbfgs), _["value"]=opt_val,
                         _["counts"]=Rcpp::NumericVector::create(1.0),
                         _["convergence"]=cc, _["message"]="lbfgs");
    } else {
      std::string alg;
      if(optimizer=="bobyqa")            alg="NLOPT_LN_BOBYQA";
      if(optimizer=="nelder-mead-nlopt") alg="NLOPT_LN_NELDERMEAD";
      if(optimizer=="subplex")           alg="NLOPT_LN_SBPLX";
      List opts = List::create(_["algorithm"]=alg, _["ftol_rel"]=reltol,
                               _["ftol_abs"]=reltol, _["xtol_rel"]=0.0001,
                               _["maxeval"]=maxit);
      List S0 = nloptr(_["x0"]=ss0, _["eval_f"]=pglmm_LL_cpp2, _["opts"]=opts,
                       _["H"]=H, _["X"]=X, _["Zt"]=Zt, _["St"]=St,
                       _["mu"]=mu, _["nested"]=nested,
                       _["REML"]=REML, _["verbose"]=verbose,
                       _["family"]=family, _["totalSize"]=totalSize);
      opt = List::create(_["par"]=S0["solution"], _["value"]=S0["objective"],
                         _["counts"]=S0["iterations"], _["convergence"]=S0["status"],
                         _["message"]=S0["message"]);
    }
    arma::vec par0 = abs(as<arma::vec>(opt["par"]));
    ss0      = wrap(par0);
    LL       = as_scalar(as<double>(opt["value"]));
    convcode = as<int>(opt["convergence"]);
    niter    = as<arma::vec>(opt["counts"]);
    est_ss   = par0; est_B = B;
    ++iteration;
    if(verbose) Rcout << "var part: " << iteration << " " << LL << std::endl;
  } // end outer PQL

  // ==== form final iV once (needed for downstream predict/simulate) ====
  if(q_Nested == 0 && q_nn > 0){
    arma::sp_mat Ut_sp = build_Ut_helper(ss0, St, Zt, q_nn);
    arma::mat Ut_d = arma::mat(Ut_sp), U_d = Ut_d.t();
    int Q = (int)Ut_d.n_rows;
    arma::vec W;
    if(family=="binomial") W = totalSize % mu % (1.0-mu);
    if(family=="poisson")  W = mu;
    arma::mat WU = U_d; WU.each_col() %= W;
    arma::mat M  = arma::inv_sympd(arma::eye(Q,Q) + Ut_d * WU);
    iV = -(WU * (M * WU.t())); iV.diag() += W;

  } else if(q_Nested > 0 && use_blk && q_nn == 0){
    // OPT6b final iV: nested-only, block-diagonal iA (no Woodbury needed)
    arma::vec pq;
    if(family=="binomial") pq = 1.0/(totalSize % mu % (1.0-mu));
    if(family=="poisson")  pq = 1.0/mu;
    arma::vec sn_vec(q_Nested);
    for(int j = 0; j < q_Nested; ++j)
      sn_vec[j] = std::abs((double)ss0[q_nn+j]);
    std::vector<arma::mat> nj_fixed(n_site,
                                    arma::mat(n_sp, n_sp, arma::fill::zeros));
    for(int j = 0; j < q_Nested; ++j){
      double snj2 = sn_vec[j] * sn_vec[j];
      sp_mat njsp = nested[j];
      for(int k = 0; k < n_site; ++k){
        int s = k*n_sp;
        nj_fixed[k] += snj2 *
          arma::mat(njsp(arma::span(s,s+n_sp-1), arma::span(s,s+n_sp-1)));
      }
    }
    std::vector<arma::mat> chols(n_site);
    block_chol(chols, pq, nj_fixed, n_sp, n_site);
    arma::mat iA_dense(n, n, arma::fill::zeros);
    for(int k = 0; k < n_site; ++k){
      int s = k*n_sp;
      arma::mat Ri = arma::inv(arma::trimatu(chols[k]));
      iA_dense(arma::span(s,s+n_sp-1), arma::span(s,s+n_sp-1)) = Ri*Ri.t();
    }
    iV = iA_dense;
  } else {
    // Fallback: q_Nested>0 with q_nn>0 (nested+non-nested), or nested-only
    // without block structure. Build A, Cholesky, dense iA via triangular
    // solves, then dense Woodbury — no explicit inv(), no sparse conversion,
    // no R-level List round-trip.
    arma::vec pq_fin;
    if(family=="binomial") pq_fin = 1.0/(totalSize % mu % (1.0-mu));
    if(family=="poisson")  pq_fin = 1.0/mu;

    arma::vec sn_fin(q_Nested);
    for(int j = 0; j < q_Nested; ++j)
      sn_fin[j] = std::abs((double)ss0[q_nn + j]);

    arma::sp_mat A_sp;
    if(family=="binomial") A_sp = sp_mat(diagmat(pq_fin));
    if(family=="poisson")  A_sp = sp_mat(diagmat(1.0/mu));
    for(int j = 0; j < q_Nested; ++j){
      double snj2 = sn_fin[j] * sn_fin[j];
      sp_mat nj = nested[j];
      A_sp = A_sp + snj2 * nj;
    }
    arma::mat A1(A_sp);
    arma::mat chol_A;
    arma::mat iA;
    if (arma::chol(chol_A, A1)) {
      arma::mat Rlo = chol_A.t();   // lower L s.t. L*L' = A
      // iA = A^{-1} via two batch triangular solves — O(n^3), no explicit inv()
      iA = arma::solve(arma::trimatl(Rlo), arma::eye(n, n));
      iA = arma::solve(arma::trimatu(chol_A), iA);
    } else {
      iA = arma::inv(A1);           // LU fallback (A not PD for final params)
    }

    if(q_nn > 0){
      arma::sp_mat Ut_sp = build_Ut_helper(ss0, St, Zt, q_nn);
      arma::mat Ut_d = arma::mat(Ut_sp), U_d = Ut_d.t();
      int Q = (int)Ut_d.n_rows;
      // Dense Woodbury: iV = iA - (iA*U) * M * (iA*U)^T,  M=(I+Ut*iA*U)^{-1}
      arma::mat B_fin = iA * U_d;                                   // n×Q
      arma::mat M_fin = arma::inv_sympd(arma::eye(Q,Q) + Ut_d * B_fin);
      iV = iA - B_fin * (M_fin * B_fin.t());
    } else {
      iV = iA;   // q_nn==0, nested-only non-block: iV = iA
    }
  }

  return List::create(_["B"]=B, _["ss"]=ss0, _["iV"]=iV,
                      _["mu"]=mu, _["H"]=H,
                      _["convcode"]=convcode, _["niter"]=niter, _["LL"]=LL);
}

// [[Rcpp::export]]
int sexp_type(SEXP x){ return TYPEOF(x); }
