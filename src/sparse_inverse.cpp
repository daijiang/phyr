// Selected inverse of a sparse symmetric matrix from its simplicial Cholesky
// factorization (CHOLMOD, as returned by Matrix::Cholesky(super = FALSE)).
//
// With P A P' = L D L' (L unit lower triangular; for an LL' factorization,
// L D^1/2), the entries of Z = (P A P')^-1 on the pattern of L + L' satisfy
// (Takahashi, Fagan and Chin 1973; Erisman and Tinney 1975)
//   Z_ij = - sum_{k in S_j} Z_ik L_kj           (i in S_j)
//   Z_jj = 1 / D_j - sum_{k in S_j} L_kj Z_kj
// with S_j the rows of the off-diagonal entries of column j of L, computed
// from the last column to the first. The pattern of L is closed for these
// sums (the rows of S_j that are larger than k are in column k), so only the
// entries of Z on the pattern of L are needed, at about the cost of the
// factorization. The pattern of A is included in that of L + L'. D can have
// negative values (LDL' factorization of a quasi-definite matrix).

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <utility>
#include <vector>
#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

// Entries of A^-1 at the positions (qi, qj) of P A P' (0-based), from the
// slots of the factor: p, i, x, nz (CHOLMOD simplicial layout: column j has
// the entries p[j] to p[j] + nz[j] - 1, with sorted rows, the diagonal first)
// and is_ll (LL' factorization: the diagonal has L_jj; otherwise D_j).
// [[Rcpp::export]]
NumericVector sparse_inv_subset(const IntegerVector& p, const IntegerVector& i,
                                const NumericVector& x, const IntegerVector& nz, bool is_ll,
                                const IntegerVector& qi, const IntegerVector& qj) {
  const int n = nz.size();
  const R_xlen_t nnz = x.size();
  std::vector<double> lx(nnz), d(n), zx(nnz), acc(n, 0.0);
  // raw pointers in the loops (element access of Rcpp vectors is checked)
  const int* pp = INTEGER(p);
  const int* pnz = INTEGER(nz);
  const int* ii = INTEGER(i);
  const double* xx = REAL(x);
  // unit L and D
  for (int j = 0; j < n; j++) {
    const int pj = pp[j];
    if (pnz[j] < 1 || ii[pj] != j) stop("sparse_inv_subset: the factor has no diagonal entry");
    const double dj = xx[pj];
    if (is_ll) {
      d[j] = dj * dj;
      for (int t = pj + 1; t < pj + pnz[j]; t++) lx[t] = xx[t] / dj;
    } else {
      d[j] = dj;
      for (int t = pj + 1; t < pj + pnz[j]; t++) lx[t] = xx[t];
    }
    if (d[j] == 0) stop("sparse_inv_subset: singular factor");
  }
  for (int j = n - 1; j >= 0; j--) {
    const int pj = pp[j], ej = pp[j] + pnz[j];
    // acc[r] = sum_{k in S_j} Z_rk L_kj for r in S_j: each pair (r, k) with
    // r >= k is read once from column k of Z, and used for Z_rk and Z_kr
    for (int t = pj + 1; t < ej; t++) {
      const int k = ii[t];
      const double lk = lx[t];
      const int ek = pp[k] + pnz[k];
      int a = pp[k];  // column k of Z (rows >= k, the diagonal first)
      // the rows of S_j from k on are in column k: merge (a linear scan is
      // faster than searching in practice)
      for (int b = t; b < ej; b++) {
        const int rb = ii[b];
        while (a < ek && ii[a] < rb) a++;
        if (a >= ek || ii[a] != rb) stop("sparse_inv_subset: the pattern of the factor is not closed");
        const double z = zx[a];
        acc[rb] += z * lk;
        if (rb > k) acc[k] += z * lx[b];
        a++;
      }
    }
    double s = 0;
    for (int t = pj + 1; t < ej; t++) {
      const int r = ii[t];
      zx[t] = -acc[r];
      acc[r] = 0;
      s += lx[t] * zx[t];
    }
    zx[pj] = 1 / d[j] - s;
  }
  const R_xlen_t nq = qi.size();
  NumericVector out(nq);
  for (R_xlen_t q = 0; q < nq; q++) {
    int r = qi[q], c = qj[q];
    if (r < c) std::swap(r, c);
    if (c < 0 || r >= n) stop("sparse_inv_subset: index out of range");
    // binary search of row r in column c
    int lo = p[c], hi = p[c] + nz[c] - 1;
    double v = NA_REAL;
    while (lo <= hi) {
      const int mid = lo + (hi - lo) / 2;
      if (i[mid] == r) {
        v = zx[mid];
        break;
      }
      if (i[mid] < r) lo = mid + 1; else hi = mid - 1;
    }
    if (ISNA(v)) stop("sparse_inv_subset: a requested entry is not in the pattern of the factor");
    out[q] = v;
  }
  return out;
}

// The same from a supernodal factorization P A P' = L L' (CHOLMOD, as returned
// by Matrix::Cholesky(super = TRUE)): slots super, pi, px, s and x. Supernode
// J has the columns c = super[J] to super[J + 1] - 1 and the rows s[pi[J]] to
// s[pi[J + 1] - 1] (the columns c first, then the other rows r, sorted), its
// values a dense column-major block at x[px[J]] (with leading dimension the
// number of rows). With Y = L_rc L_cc^-1, the recursion of Takahashi for the
// columns of J is, with dense blocks (Lin et al. 2011, SelInv):
//   Z_rc = - Z_rr Y
//   Z_cc = (L_cc L_cc')^-1 - Y' Z_rc
// from the last supernode to the first. Z_rr is on the pattern of the
// ancestors of J (the pattern is closed, as for the simplicial factor: it is
// what the supernodal factorization uses too). Z is stored as L (the
// explicit zeros of relaxed supernodes only make the blocks larger).
// [[Rcpp::export]]
NumericVector sparse_inv_subset_super(const IntegerVector& super, const IntegerVector& pi,
                                      const IntegerVector& px, const IntegerVector& s,
                                      const NumericVector& x, const IntegerVector& qi,
                                      const IntegerVector& qj) {
  const int ns = super.size() - 1;
  const int n = super[ns];
  const int* sup = INTEGER(super);
  const int* ppi = INTEGER(pi);
  const int* ppx = INTEGER(px);
  const int* ss = INTEGER(s);
  const double* xx = REAL(x);
  std::vector<int> col2s(n), relpos(n, -1);
  int mmax = 0, ymax = 0;
  for (int J = 0; J < ns; J++) {
    const int nc = sup[J + 1] - sup[J], m = ppi[J + 1] - ppi[J] - nc;
    for (int c = sup[J]; c < sup[J + 1]; c++) col2s[c] = J;
    if (m > mmax) mmax = m;
    if (m * nc > ymax) ymax = m * nc;
  }
  std::vector<double> z(x.size()), zrr((size_t) mmax * mmax), y(ymax);
  const double one = 1, mone = -1, zero = 0;
  for (int J = ns - 1; J >= 0; J--) {
    const int nc = sup[J + 1] - sup[J], nr = ppi[J + 1] - ppi[J], m = nr - nc;
    const double* lx = xx + ppx[J];
    double* zx = z.data() + ppx[J];
    // (L_cc L_cc')^-1 in the lower triangle of the diagonal block of Z
    for (int j = 0; j < nc; j++) {
      for (int t = j; t < nc; t++) zx[t + (size_t) j * nr] = lx[t + (size_t) j * nr];
    }
    int info = 0;
    F77_CALL(dpotri)("L", &nc, zx, &nr, &info FCONE);
    if (info != 0) stop("sparse_inv_subset_super: singular factor");
    if (m == 0) continue;
    // Y = L_rc L_cc^-1 (m by nc)
    for (int j = 0; j < nc; j++) {
      for (int t = 0; t < m; t++) y[t + (size_t) j * m] = lx[nc + t + (size_t) j * nr];
    }
    F77_CALL(dtrsm)("R", "L", "N", "N", &m, &nc, &one, lx, &nr, y.data(), &m
                    FCONE FCONE FCONE FCONE);
    // Z_rr (lower triangle) from the columns of the ancestors: the rows r
    // from row r_a on are in the column of r_a
    const int* rows = ss + ppi[J] + nc;
    int K = -1;
    for (int a = 0; a < m; a++) {
      const int ra = rows[a];
      if (col2s[ra] != K) {
        K = col2s[ra];
        for (int t = ppi[K]; t < ppi[K + 1]; t++) relpos[ss[t]] = t - ppi[K];
      }
      const int nrK = ppi[K + 1] - ppi[K];
      const double* zk = z.data() + ppx[K] + (size_t) (ra - sup[K]) * nrK;
      for (int b = a; b < m; b++) {
        const int pos = relpos[rows[b]];
        if (pos < 0 || pos >= nrK || ss[ppi[K] + pos] != rows[b]) {
          stop("sparse_inv_subset_super: the pattern of the factor is not closed");
        }
        zrr[b + (size_t) a * m] = zk[pos];
      }
    }
    // Z_rc = - Z_rr Y, then Z_cc -= Y' Z_rc
    F77_CALL(dsymm)("L", "L", &m, &nc, &mone, zrr.data(), &m, y.data(), &m, &zero, zx + nc, &nr
                    FCONE FCONE);
    F77_CALL(dgemm)("T", "N", &nc, &nc, &m, &mone, y.data(), &m, zx + nc, &nr, &one, zx, &nr
                    FCONE FCONE);
  }
  const R_xlen_t nq = qi.size();
  NumericVector out(nq);
  for (R_xlen_t q = 0; q < nq; q++) {
    int r = qi[q], c = qj[q];
    if (r < c) std::swap(r, c);
    if (c < 0 || r >= n) stop("sparse_inv_subset_super: index out of range");
    const int K = col2s[c], jk = c - sup[K], nrK = ppi[K + 1] - ppi[K];
    // binary search of row r in column c (rows from the diagonal on)
    int lo = jk, hi = nrK - 1, pos = -1;
    while (lo <= hi) {
      const int mid = lo + (hi - lo) / 2;
      const int rm = ss[ppi[K] + mid];
      if (rm == r) {
        pos = mid;
        break;
      }
      if (rm < r) lo = mid + 1; else hi = mid - 1;
    }
    if (pos < 0) stop("sparse_inv_subset_super: a requested entry is not in the pattern of the factor");
    out[q] = z[ppx[K] + (size_t) jk * nrK + pos];
  }
  return out;
}
