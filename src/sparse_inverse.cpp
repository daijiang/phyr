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

#include <Rcpp.h>
#include <cmath>
#include <utility>
#include <vector>

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
