// pglmm_utils.h — shared helpers for block-diagonal Cholesky used by both
// pglmm_binary.cpp (GLMM inner loop) and pglmm_gaussian.cpp (LMM LL function).
//
// All functions are inline to avoid ODR violations when included in multiple TUs.
//
// [[Rcpp::depends(RcppArmadillo)]]
#pragma once
#include "RcppArmadillo.h"

// Count nnz in row 0 of a sparse matrix = block size for kron(I_ns, Vphy).
// Uses element access (avoids const_row_iterator portability issues).
inline int detect_n_sp(const arma::sp_mat& nj) {
  int cnt = 0;
  for (arma::uword c = 0; c < nj.n_cols; ++c) {
    if (nj.at(0, c) != 0.0) ++cnt;
    else if (cnt > 0) break;   // first zero after non-zeros = end of block
  }
  return cnt;
}

// Verify that nj is truly block-diagonal with blocks of size n_sp.
// Two fast checks (both O(nnz) or better):
//   1. Total nnz == n_site * n_sp^2  (rules out interleaved / ragged patterns)
//   2. Column n_sp (first col of block 1) has no non-zeros in rows [0, n_sp)
//      (column iterators are sorted by row, so this is a single comparison)
// Returns false → fall back to full dense path; true → safe to use block Chol.
inline bool verify_block_structure(const arma::sp_mat& nj, int n_sp) {
  arma::uword n         = nj.n_rows;
  arma::uword nsp_u     = (arma::uword)n_sp;
  arma::uword n_site_u  = n / nsp_u;

  // Check 1: exact nnz count for a pure block-diagonal matrix
  if (nj.n_nonzero != n_site_u * nsp_u * nsp_u) return false;

  // Trivially block-diagonal when there is only one block
  if (n_site_u < 2) return true;

  // Check 2: first non-zero row in column n_sp must be >= n_sp
  arma::sp_mat::const_col_iterator it = nj.begin_col(nsp_u);
  if (it != nj.end_col(nsp_u) && (int)it.row() < n_sp) return false;

  return true;
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
    if (!arma::chol(chols[k], Ak))
      arma::chol(chols[k], Ak + 1e-10 * arma::eye(n_sp, n_sp));
    logdetA += 2.0 * arma::accu(arma::log(chols[k].diag()));
  }
  return logdetA;
}
