// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-
// Binomial / poisson PGLMM fitted by penalized quasi-likelihood (PQL), for
// models whose nested terms are absent or all DIAGONAL (e.g. the
// observation-level term "1|obs" that is added by default). The algorithm and
// the structure of the PQL loop are the same as in pglmm_binary.cpp, but:
//  (1) diagonal nested terms are folded into the diagonal part of V, so the
//      Woodbury identity is used instead of a dense n x n Cholesky;
//  (2) every quantity is formed with Zu, the unique columns of Zt
//      (observations of the same levels have the same column):
//      Zt diag(w) Zt' = Zu diag(groupsum(w)) Zu'; there is no dense Q x n or
//      n x n work;
//  (3) when the weights do not depend on the parameters (no diagonal nested
//      term), the statistics are built once per optimizer call and each
//      likelihood evaluation is one Q x Q Cholesky;
//  (4) the derivative-free optimizer is called with the NLopt C API.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
#include "RcppArmadillo.h"
#include <nloptrAPI.h>
#include <unordered_map>
#include <string>
#include <vector>
using namespace Rcpp;

struct QB {
  arma::mat Zu;            // Q x n_u : unique columns of Zt
  arma::uvec grp;          // n : column of Zu for each observation
  arma::mat St_d;          // q_nn x Q
  std::vector<arma::vec> dn;   // diagonals of the (diagonal) nested terms
  int n, p, Q, q_nn, q_dn, n_u;
  bool REML;
  // per-optimiser-call data
  const arma::mat* X; const arma::vec* H; arma::vec iW;   // iW = 1/W
  bool have_fixed; arma::mat G_fixed; arma::vec ZWH_fixed; arma::mat ZWX_fixed, XWX_fixed; double HWH_fixed, sumloga_fixed;
};

static arma::vec gsum(const arma::vec& x, const arma::uvec& grp, int n_u) {
  arma::vec o(n_u, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_elem; ++i) o[grp[i]] += x[i];
  return o;
}
static arma::mat gsum_m(const arma::mat& x, const arma::uvec& grp, int n_u) {
  arma::mat o(n_u, x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j)
    for (arma::uword i = 0; i < x.n_rows; ++i) o(grp[i], j) += x(i, j);
  return o;
}
static arma::mat G_of(const QB& P, const arma::vec& wu) {      // Zu diag(wu) Zu'
  arma::mat ZuW = P.Zu; ZuW.each_row() %= wu.t();
  return ZuW * P.Zu.t();
}

static void qb_setup(QB& P, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, int p, bool REML) {
  P.n = Zt.n_cols; P.Q = Zt.n_rows; P.q_nn = St.n_rows; P.p = p; P.REML = REML;
  P.St_d = arma::mat(St);
  // unique columns of Zt
  std::unordered_map<std::string, arma::uword> seen;
  P.grp.set_size(P.n);
  std::vector<arma::uword> first;
  for (int j = 0; j < P.n; ++j) {
    std::string key;
    for (arma::sp_mat::const_col_iterator it = Zt.begin_col(j); it != Zt.end_col(j); ++it) {
      arma::uword r = it.row(); double v = *it;
      key.append((const char*)&r, sizeof(r)); key.append((const char*)&v, sizeof(v));
    }
    std::unordered_map<std::string, arma::uword>::iterator f = seen.find(key);
    if (f == seen.end()) { arma::uword id = first.size(); seen.emplace(key, id); first.push_back(j); P.grp[j] = id; }
    else P.grp[j] = f->second;
  }
  P.n_u = (int)first.size();
  P.Zu.zeros(P.Q, P.n_u);
  for (int g = 0; g < P.n_u; ++g)
    for (arma::sp_mat::const_col_iterator it = Zt.begin_col(first[g]); it != Zt.end_col(first[g]); ++it) P.Zu(it.row(), g) = *it;
  // diagonal nested terms
  P.q_dn = nested.size();
  for (int j = 0; j < P.q_dn; ++j) {
    arma::sp_mat nj = nested[j];
    arma::vec d(nj.diag());
    if ((int)nj.n_nonzero != (int)arma::accu(d != 0.0)) Rcpp::stop("nested terms must be diagonal here");
    P.dn.push_back(d);
  }
  P.have_fixed = false;
}

// LL of the PQL variance step (same quantity as pglmm_LL_cpp).
static double qb_ll(const double* x, QB& P) {
  arma::rowvec sr(P.q_nn);
  for (int r = 0; r < P.q_nn; ++r) sr[r] = std::abs(x[r]);
  arma::vec c = arma::vectorise(sr * P.St_d);
  const arma::vec& H = *P.H; const arma::mat& X = *P.X;

  arma::mat G, ZWX, XWX; arma::vec ZWH; double HWH, sumloga;
  if (P.q_dn == 0 && P.have_fixed) {
    G = P.G_fixed; ZWH = P.ZWH_fixed; HWH = P.HWH_fixed; sumloga = P.sumloga_fixed;
    if (P.REML) { ZWX = P.ZWX_fixed; XWX = P.XWX_fixed; }
  } else {
    arma::vec a = P.iW;
    for (int j = 0; j < P.q_dn; ++j) { double s = x[P.q_nn + j]; a += (s * s) * P.dn[j]; }
    arma::vec w = 1.0 / a;
    G   = G_of(P, gsum(w, P.grp, P.n_u));
    ZWH = P.Zu * gsum(w % H, P.grp, P.n_u);
    HWH = arma::dot(w % H, H);
    sumloga = arma::accu(arma::log(a));
    if (P.REML) {
      arma::mat WX = X; WX.each_col() %= w;
      ZWX = P.Zu * gsum_m(WX, P.grp, P.n_u);
      XWX = X.t() * WX;
    }
    if (P.q_dn == 0) {   // cache: weights are parameter-free
      P.G_fixed = G; P.ZWH_fixed = ZWH; P.HWH_fixed = HWH; P.sumloga_fixed = sumloga;
      if (P.REML) { P.ZWX_fixed = ZWX; P.XWX_fixed = XWX; }
      P.have_fixed = true;
    }
  }
  arma::mat K = G; K.each_col() %= c; K.each_row() %= c.t(); K.diag() += 1.0;
  arma::mat R;
  if (!arma::chol(R, arma::symmatu(K))) return 1e15;
  double logdetV = 2.0 * arma::accu(arma::log(R.diag())) + sumloga;
  arma::vec u = arma::solve(arma::trimatl(R.t()), c % ZWH);       // R' u = Ut W H
  double HiVH = HWH - arma::dot(u, u);
  double LL = 0.5 * (logdetV + HiVH);
  if (P.REML) {
    arma::mat cZWX = ZWX; cZWX.each_col() %= c;
    arma::mat T = arma::solve(arma::trimatl(R.t()), cZWX);
    double ld, sg; arma::log_det(ld, sg, XWX - T.t() * T);
    LL += 0.5 * ld;
  }
  return LL;
}
static double qb_cb(unsigned npar, const double* x, double* grad, void* data) {
  QB* P = static_cast<QB*>(data);
  double v;
  try { v = qb_ll(x, *P); } catch (...) { v = 1e30; }
  return v;
}

// [[Rcpp::export]]
List pglmm_internal_q_cpp(const arma::mat& X, const arma::vec& Y, const arma::sp_mat& Zt, const arma::sp_mat& St,
                 const List& nested, const bool REML, const int maxit, const double reltol,
                 const double tol_pql, const double maxit_pql, const std::string optimizer,
                 arma::mat B_init, arma::vec ss, const std::string family, arma::vec totalSize) {
  int n = X.n_rows, p = X.n_cols, q = ss.n_elem;
  QB P; qb_setup(P, Zt, St, nested, p, REML);
  int q_nn = P.q_nn;
  arma::vec B = B_init.col(0), b(n, arma::fill::zeros), mu, Z, H, W;
  auto linkinv = [&](const arma::vec& eta) -> arma::vec {
    if (family == "binomial") return arma::vec(arma::exp(eta) / (1 + arma::exp(eta)));
    return arma::vec(arma::exp(eta));
  };
  auto workW = [&]() -> arma::vec { return (family == "binomial") ? arma::vec(totalSize % mu % (1.0 - mu)) : mu; };
  auto workZ = [&]() -> arma::vec {
    if (family == "binomial") return arma::vec(X * B + b + (Y / totalSize - mu) / (mu % (1.0 - mu)));
    return arma::vec(X * B + b + (Y - mu) / mu);
  };
  mu = linkinv(X * B);
  arma::vec est_ss = ss, est_B = B, oldest_ss(q), oldest_B(p); oldest_ss.fill(1e6); oldest_B.fill(1e6);
  std::vector<double> ss0(ss.begin(), ss.end());
  double tol2 = tol_pql * tol_pql, LL = 0.0;
  int iteration = 0, convcode = 0;
  arma::vec c, a, w; arma::mat Kinv;

  // state for the current (ss0, mu): c, a, w, Kinv
  auto refresh = [&]() {
    arma::rowvec sr(q_nn); for (int r = 0; r < q_nn; ++r) sr[r] = ss0[r];
    c = arma::vectorise(sr * P.St_d);
    a = 1.0 / workW();
    for (int j = 0; j < P.q_dn; ++j) a += (ss0[q_nn + j] * ss0[q_nn + j]) * P.dn[j];
    w = 1.0 / a;
    arma::mat K = G_of(P, gsum(w, P.grp, P.n_u)); K.each_col() %= c; K.each_row() %= c.t(); K.diag() += 1.0;
    Kinv = arma::inv_sympd(arma::symmatu(K));
  };
  auto Ut_mult = [&](const arma::mat& x) -> arma::mat {        // Ut x  (Q x k)
    arma::mat o = P.Zu * gsum_m(x, P.grp, P.n_u); o.each_col() %= c; return o;
  };
  auto U_mult = [&](const arma::vec& y) -> arma::vec {         // U y   (n)
    arma::vec t = P.Zu.t() * (c % y); arma::vec o(n);
    for (int i = 0; i < n; ++i) o[i] = t[P.grp[i]];
    return o;
  };
  auto iV_mult = [&](const arma::mat& x) -> arma::mat {        // iV x  (n x k)
    arma::mat wx = x; wx.each_col() %= w;
    arma::mat t = Kinv * Ut_mult(wx);
    arma::mat o = wx;
    for (arma::uword k = 0; k < x.n_cols; ++k) o.col(k) -= w % U_mult(t.col(k));
    return o;
  };

  while ((arma::dot(est_ss - oldest_ss, est_ss - oldest_ss) > tol2 || arma::dot(est_B - oldest_B, est_B - oldest_B) > tol2) &&
         iteration <= maxit_pql) {
    oldest_ss = est_ss; oldest_B = est_B;
    Rcpp::checkUserInterrupt();
    int iteration_m = 0; arma::vec est_B_m = B, oldest_B_m(p); oldest_B_m.fill(1e6);
    while (arma::dot(est_B_m - oldest_B_m, est_B_m - oldest_B_m) > tol2 && iteration_m <= maxit_pql) {
      oldest_B_m = est_B_m;
      refresh();
      Z = workZ();
      arma::mat iVX = iV_mult(X);
      arma::vec iVZ = iV_mult(arma::mat(Z)).col(0);
      B = arma::solve(X.t() * iVX, X.t() * iVZ);
      arma::vec r = Z - X * B;
      arma::vec iVr = iV_mult(arma::mat(r)).col(0);
      b = U_mult(Ut_mult(arma::mat(iVr)).col(0));
      for (int j = 0; j < P.q_dn; ++j) b += (ss0[q_nn + j] * ss0[q_nn + j]) * (P.dn[j] % iVr);
      mu = linkinv(X * B + b);
      est_B_m = B; ++iteration_m;
      if (B.has_nan()) Rcpp::stop("Estimation of B failed.");
    }

    Z = workZ(); H = Z - X * B;
    P.X = &X; P.H = &H; P.iW = 1.0 / workW(); P.have_fixed = false;
    nlopt_algorithm alg = NLOPT_LN_NELDERMEAD;
    if (optimizer == "bobyqa") alg = NLOPT_LN_BOBYQA;
    if (optimizer == "subplex") alg = NLOPT_LN_SBPLX;
    nlopt_opt opt = nlopt_create(alg, (unsigned)q);
    nlopt_set_min_objective(opt, qb_cb, &P);
    nlopt_set_ftol_rel(opt, reltol); nlopt_set_ftol_abs(opt, reltol);
    nlopt_set_xtol_rel(opt, 0.0001); nlopt_set_maxeval(opt, maxit);
    double val = 1e30;
    nlopt_result res = nlopt_optimize(opt, ss0.data(), &val);
    nlopt_destroy(opt);
    convcode = (int)res; LL = val;
    for (int k = 0; k < q; ++k) ss0[k] = std::abs(ss0[k]);
    est_ss = arma::vec(ss0); est_B = B; ++iteration;
  }

  refresh();                        // factors at the final (ss, mu)
  return List::create(_["B"] = arma::mat(B), _["ss"] = arma::mat(arma::vec(ss0)), _["mu"] = arma::mat(mu),
                      _["H"] = H, _["convcode"] = convcode, _["LL"] = LL,
                      _["niter"] = arma::mat(1, 1, arma::fill::ones),
                      _["B.cov"] = arma::inv_sympd(arma::symmatu(X.t() * iV_mult(X))));
}
