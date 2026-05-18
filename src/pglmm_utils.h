// pglmm_utils.h — shared helpers for block-diagonal Cholesky used by both
// pglmm_binary.cpp (GLMM inner loop) and pglmm_gaussian.cpp (LMM LL function).
//
// All functions are inline to avoid ODR violations when included in multiple TUs.
//
// [[Rcpp::depends(RcppArmadillo)]]
#pragma once
#include "RcppArmadillo.h"

// Count consecutive NZs in row 0 of a sparse matrix.
// May under-count when Vphy has internal zero off-diagonals (species diverging
// at root give exact zeros in vcv, making row 0 non-contiguous).
// Use find_block_size() for a robust detection that handles this case.
inline int detect_n_sp(const arma::sp_mat& nj) {
  int cnt = 0;
  for (arma::uword c = 0; c < nj.n_cols; ++c) {
    if (nj.at(0, c) != 0.0) ++cnt;
    else if (cnt > 0) break;
  }
  return cnt;
}

// Verify that nj is truly block-diagonal with blocks of size n_sp.
// Three O(1) checks on boundary columns; does NOT require exactly
// n_site * n_sp^2 non-zeros (Vphy can have zero off-diagonals for species
// that diverge at the root, so block columns may be sparse).
// Returns false → fall back to full dense path; true → safe to use block Chol.
inline bool verify_block_structure(const arma::sp_mat& nj, int n_sp) {
  arma::uword n        = nj.n_rows;
  arma::uword nsp_u    = (arma::uword)n_sp;
  arma::uword n_site_u = n / nsp_u;

  // Trivially block-diagonal when there is only one block.
  if (n_site_u < 2) return true;

  // Check A: first NZ row in column n_sp (first col of block 1) must be >= n_sp.
  {
    arma::sp_mat::const_col_iterator it = nj.begin_col(nsp_u);
    if (it != nj.end_col(nsp_u) && (int)it.row() < n_sp) return false;
  }

  // Check B: all NZs in column n_sp-1 (last col of block 0) must have row < n_sp.
  {
    arma::sp_mat::const_col_iterator it = nj.begin_col(nsp_u - 1);
    while (it != nj.end_col(nsp_u - 1)) {
      if ((int)it.row() >= n_sp) return false;
      ++it;
    }
  }

  // Check C: first NZ in last column must be >= n - n_sp.
  // For kron(I_site, Vphy), last col belongs to last block (rows [n-n_sp, n)).
  // For a wrong n_sp_cand < true_n_sp: first_nz_row = n - true_n_sp < n - n_sp_cand,
  // so this check reliably rejects under-estimated block sizes.
  {
    arma::sp_mat::const_col_iterator it = nj.begin_col(n - 1);
    if (it != nj.end_col(n - 1) && (int)it.row() < (int)(n - nsp_u)) return false;
  }

  return true;
}

// Robustly find the Kronecker block size for nj = kron(I_n_site, Vphy).
// Strategy: start with detect_n_sp() as a cheap candidate; if it fails
// verify_block_structure(), search upward through divisors of n.
// verify_block_structure()'s Check C guarantees that any k < true_n_sp is
// rejected (first NZ in last col = n - true_n_sp < n - k), so the search
// stops at exactly the true block size.
// Returns 0 if no valid block size is found (fall back to full dense path).
inline int find_block_size(const arma::sp_mat& nj, int n) {
  int cand = detect_n_sp(nj);
  if (cand >= 2 && n % cand == 0 && verify_block_structure(nj, cand))
    return cand;

  // Search upward: for typical community data (n_sp ~ 5-100),
  // this loop runs at most n_sp iterations, each doing 3 column lookups.
  int start = std::max(cand + 1, 2);
  for (int k = start; k <= n / 2; ++k) {
    if (n % k == 0 && verify_block_structure(nj, k))
      return k;
  }
  return 0;
}

// Compute A^{-1} from upper Cholesky factor R where A = R'R.
// Polyfill for arma::chol2inv (not available in all Armadillo versions).
// Uses back-substitution: iR = R^{-1}, then A^{-1} = iR * iR'.
inline arma::mat phyr_chol2inv(const arma::mat& R) {
  int n = (int)R.n_rows;
  arma::mat iR = arma::solve(arma::trimatu(R), arma::eye<arma::mat>(n, n));
  return arma::symmatu(iR * iR.t());
}

// Compute block Cholesky factors for A_k = diag(pq_k) + nj_fixed[k];
// returns logdetA = sum_k 2*sum(log(diag(R_k))).
// On entry:  chols is a pre-allocated vector<arma::mat> of length n_site.
// On return: chols[k] holds the upper Cholesky factor R_k for block k.
inline double block_chol(std::vector<arma::mat>& chols,
                         const arma::vec& pq,
                         const std::vector<arma::mat>& nj_fixed,
                         int n_sp, int n_site) {
  double logdetA = 0.0;
  for (int k = 0; k < n_site; ++k) {
    int s = k * n_sp;
    arma::mat Ak = arma::diagmat(pq.subvec(s, s + n_sp - 1)) + nj_fixed[k];
    if (!arma::chol(chols[k], Ak)) {
      arma::mat Ak_j = Ak + 1e-10 * arma::eye(n_sp, n_sp);
      if (!arma::chol(chols[k], Ak_j))
        arma::chol(chols[k], Ak_j + 1e-6 * arma::eye(n_sp, n_sp));
    }
    logdetA += 2.0 * arma::accu(arma::log(chols[k].diag()));
  }
  return logdetA;
}
