// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef __PHYR_COR_PHYLO_H
#define __PHYR_COR_PHYLO_H

#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>
#include <vector>


typedef uint_fast32_t uint_t;

using namespace Rcpp;


/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Classes
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */


// Matrices to be kept in output
class cp_matrices {
public:
  arma::mat mean_sd_X;
  std::vector<arma::vec> sd_U;
  arma::mat XX;
  arma::mat UU;
  arma::mat MM;
  arma::mat Vphy;
  arma::mat R;
  arma::mat V;
  arma::mat C;
  arma::mat iD;
  arma::mat B;
  
  cp_matrices() {}
  cp_matrices(const arma::mat& mean_sd_X_, const std::vector<arma::vec>& sd_U_, 
              const arma::mat& XX_, const arma::mat& UU_, const arma::mat& MM_, 
              const arma::mat& Vphy_, const arma::mat& R_, const arma::mat& V_, 
              const arma::mat& C_, const arma::mat& B_) 
    : mean_sd_X(mean_sd_X_), sd_U(sd_U_), XX(XX_), UU(UU_), MM(MM_), Vphy(Vphy_), 
      R(R_), V(V_), C(C_), B(B_) {
    iD = arma::chol(V).t();
  };
  
  
};



// Info to calculate the log-likelihood
class LL_info {
public:
  arma::vec par0;  // par to start with
  arma::mat XX;
  arma::mat UU;
  arma::mat MM;
  arma::mat Vphy;
  arma::mat tau;
  bool REML;
  bool constrain_d;
  bool verbose;
  uint_t iters;
  arma::vec min_par; // par for minimum LL
  double LL;
  int convcode;
  
  LL_info() {}
  LL_info(const arma::mat& X,
          const std::vector<arma::mat>& U,
          const arma::mat& M,
          const arma::mat& Vphy_,
          const bool& REML_,
          const bool& constrain_d_,
          const bool& verbose_);
  
  // LL_info(cp_matrices cpm, const bool& REML_, const bool& constrain_d_, 
  //        const bool& verbose_);
  
};


class boot_results {
public:
  arma::cube corrs;
  arma::mat d;
  arma::mat B0;
  arma::cube B_cov;
  arma::cube X_means_sds;
  
  boot_results(const uint_t& p, const uint_t& n_covs, const uint_t& n_reps) 
    : corrs(p,p, n_reps), 
      d(p, n_reps), 
      B0(p + n_covs, n_reps), 
      B_cov(p + n_covs, p + n_covs, n_reps),
      X_means_sds(p, 2, n_reps) {};
  boot_results(const XPtr<cp_matrices> cpm, const uint_t& n_reps) 
    : corrs(cpm->mean_sd_X.n_rows, cpm->mean_sd_X.n_rows, n_reps), 
      d(cpm->mean_sd_X.n_rows, n_reps),
      B0(cpm->UU.n_cols, n_reps),
      B_cov(cpm->UU.n_cols, cpm->UU.n_cols, n_reps),
      X_means_sds(cpm->mean_sd_X.n_rows, 2, n_reps) {};
};





/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 Miscellaneous helper functions
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */

// Flexible power function needed for multiple functions below
// a^b

arma::vec flex_pow(const arma::vec& a, const double& b) {
  uint_t n = a.n_elem;
  arma::vec x(n);
  for (uint_t i = 0; i < n; i++) x(i) = std::pow(a(i), b);
  return x;
}
arma::mat flex_pow(const double& a, const arma::mat& b) {
  uint_t nr = b.n_rows, nc = b.n_cols;
  arma::mat x(nr, nc);
  for (uint_t i = 0; i < nr; i++) {
    for (uint_t j = 0; j < nc; j++) {
      x(i, j) = std::pow(a, b(i, j));
    }
  }
  return x;
}
arma::mat flex_pow(const arma::mat& a, const double& b) {
  uint_t nr = a.n_rows, nc = a.n_cols;
  arma::mat x(nr, nc);
  for (uint_t i = 0; i < nr; i++) {
    for (uint_t j = 0; j < nc; j++) {
      x(i, j) = std::pow(a(i, j), b);
    }
  }
  return x;
}

// Transpose functions that return what I want them to:
inline arma::mat tp(const arma::mat& M){
  return M.t();
}
inline arma::rowvec tp(const arma::vec& V){
  return arma::conv_to<arma::rowvec>::from(V.t());
}
inline arma::vec tp(const arma::rowvec& V){
  return arma::conv_to<arma::vec>::from(V.t());
}

// pnorm for standard normal (i.e., ~ N(0,1))
arma::vec pnorm_cpp(const arma::vec& values, const bool& lower_tail) {
  arma::vec out(values.n_elem);
  for (uint_t i = 0; i < values.size(); i++) {
    out(i) = 0.5 * std::erfc(-values(i) * M_SQRT1_2);
  }
  if (!lower_tail) out = 1 - out;
  return out;
}





/*
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 
 "make" functions, that create objects used in multiple functions inside corphylo.cpp
 
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 ***************************************************************************************
 */



inline arma::mat make_L(const arma::vec& par, const uint_t& n, const uint_t& p) {
  arma::vec L_elements = par(arma::span(0, (p + p * (p - 1)/2) - 1));
  arma::mat L(p, p, arma::fill::zeros);
  for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
    L(arma::span(i, p-1), i) = L_elements(arma::span(j, k));
    j = k + 1;
    k += (p - i - 1);
  }
  return L;
}

inline arma::vec make_d(const arma::vec& par, const uint_t& n, const uint_t& p,
                        const bool& constrain_d, bool do_checks) {
  arma::vec d;
  if (constrain_d) {
    arma::vec logit_d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
    if (do_checks) {
      if (arma::max(arma::abs(logit_d)) > 10) return d;
    }
    d = 1/(1 + arma::exp(-logit_d));
  } else {
    d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
    if (do_checks) {
      if (max(d) > 10) d.reset();
    }
  }
  return d;
}
inline arma::vec make_d(const arma::vec& par, const uint_t& n, const uint_t& p,
                        const bool& constrain_d) {
  arma::vec d;
  if (constrain_d) {
    arma::vec logit_d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
    d = 1/(1 + arma::exp(-logit_d));
  } else {
    d = par(arma::span((p + p * (p - 1) / 2), par.n_elem - 1));
  }
  return d;
}


// OU transform
inline arma::mat make_C(const uint_t& n, const uint_t& p,
                        const arma::mat& tau, const arma::vec& d, 
                        const arma::mat& Vphy, const arma::mat& R) {
  arma::mat C(p * n, p * n, arma::fill::zeros);
  for (uint_t i = 0; i < p; i++) {
    arma::mat Cd;
    for (uint_t j = 0; j < p; j++) {
      Cd = flex_pow(d(i), tau) % flex_pow(d(j), tp(tau)) % 
        (1 - flex_pow(d(i) * d(j), Vphy));
      Cd /= (1 - d(i) * d(j));
      C(arma::span(n * i, (i + 1) * n - 1), arma::span(n * j, (j + 1) * n - 1)) =
        R(i, j) * Cd;
    }
  }
  return C;
}

inline arma::mat make_V(const arma::mat& C, const arma::mat& MM) {
  arma::mat V = C;
  V += arma::diagmat(arma::vectorise(MM));
  return V;
}

// Correlation matrix
inline arma::mat make_corrs(const arma::mat& R) {
  arma::mat Rd = arma::diagmat(flex_pow(static_cast<arma::vec>(arma::diagvec(R)), 
                                        -0.5));
  arma::mat corrs = Rd * R * Rd;
  return corrs;
}


/*
 Make vector of standard deviations and correct matrix of coefficient estimates
 and standard errors.
 */
inline void make_sd_B_mat_cov(arma::mat& B, arma::vec& sd_vec,
                              arma::mat& B_cov, arma::vec& B0,
                              const uint_t& p, 
                              const arma::mat& iV,
                              const arma::mat& UU,
                              const arma::mat& mean_sd_X,
                              const std::vector<arma::vec>& sd_U) {
  
  sd_vec.zeros(UU.n_cols);
  for (uint_t counter = 0, i = 0; i < p; counter++, i++) {
    B0[counter] += mean_sd_X(i,0);
    sd_vec[counter] = mean_sd_X(i,1);
    if (sd_U[i].n_elem > 0) {
      const arma::vec& sd_Ui(sd_U[i]);
      for (uint_t j = 0; j < sd_Ui.n_elem; j++) {
        if (sd_Ui(j) > 0) {
          counter++;
          double sd_ratio = mean_sd_X(i,1) / sd_Ui(j);
          B0[counter] *= sd_ratio;
          sd_vec[counter] = sd_ratio;
        }
      }
    }
  }
  
  B_cov = arma::inv(tp(UU) * iV * UU);
  B_cov = arma::diagmat(sd_vec) * B_cov * arma::diagmat(sd_vec);
  
  B.set_size(B0.n_elem, 4);
  B.col(0) = B0;  // Estimates
  B.col(1) = flex_pow(static_cast<arma::vec>(arma::diagvec(B_cov)), 0.5);  // SE
  B.col(2) = B0 / B.col(1); // Z-score
  B.col(3) = 2 * pnorm_cpp(arma::abs(B.col(2)), false);  // P-value
  
  return;
}





#endif
