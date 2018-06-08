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
  // Used in bootstrapping
  LL_info(const arma::mat& X,
          const std::vector<arma::mat>& U,
          const arma::mat& M,
          const LL_info& other);
  // Copy constructor
  LL_info(const LL_info& ll_info2) {
    par0 = ll_info2.par0;
    XX = ll_info2.XX;
    UU = ll_info2.UU;
    MM = ll_info2.MM;
    Vphy = ll_info2.Vphy;
    tau = ll_info2.tau;
    REML = ll_info2.REML;
    constrain_d = ll_info2.constrain_d;
    verbose = ll_info2.verbose;
    iters = ll_info2.iters;
    min_par = ll_info2.min_par;
    LL = ll_info2.LL;
    convcode = ll_info2.convcode;
  }
  
};



// Results from bootstrapping

class boot_results {
public:
  arma::cube corrs;
  arma::mat B0;
  arma::cube B_cov;
  arma::mat d;
  std::vector<arma::mat> out_mats;
  std::vector<uint_t> out_inds;
  std::vector<int> out_codes;

  boot_results(const uint_t& p, const uint_t& B_rows, const uint_t& n_reps) 
    : corrs(p, p, n_reps, arma::fill::zeros), 
      B0(B_rows, n_reps, arma::fill::zeros), 
      B_cov(B_rows, B_rows, n_reps, arma::fill::zeros),
      d(p, n_reps, arma::fill::zeros), 
      out_mats(), out_inds(), out_codes() {};

  // Insert values into a boot_results object
  void insert_values(const uint_t& i,
                     const arma::mat& corrs_i,
                     const arma::vec& B0_i,
                     const arma::mat& B_cov_i,
                     const arma::vec& d_i) {
    
    corrs.slice(i) = corrs_i;
    B0.col(i) = B0_i;
    B_cov.slice(i) = B_cov_i;
    d.col(i) = d_i;
    return;
    
  }
  
};


/*
 Matrices to be kept for bootstrapping
 One per core if doing multi-threaded
 */
class boot_mats {
public:
  // original input matrices
  const arma::mat X;  
  const std::vector<arma::mat> U;
  const arma::mat M;
  arma::mat X_new;
  
  boot_mats(const arma::mat& X_, const std::vector<arma::mat>& U_,
            const arma::mat& M_,
            const arma::mat& B_, const arma::vec& d_, const LL_info& ll_info);
  
  XPtr<LL_info> iterate(LL_info& ll_info);
  
  void one_boot(XPtr<LL_info>& ll_info_xptr, boot_results& br,
                const uint_t& i, const double& rel_tol, const int& max_iter,
                const std::string& method, const std::string& keep_boots,
                const std::vector<double>& sann);
  
  
private:
  arma::mat iD;
  arma::mat X_pred;

  // Method for returning bootstrapped data
  void boot_data(LL_info& ll_info, boot_results& br, const uint_t& i);

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

inline arma::vec flex_pow(const arma::vec& a, const double& b) {
  uint_t n = a.n_elem;
  arma::vec x(n);
  for (uint_t i = 0; i < n; i++) x(i) = std::pow(a(i), b);
  return x;
}
inline arma::mat flex_pow(const double& a, const arma::mat& b) {
  uint_t nr = b.n_rows, nc = b.n_cols;
  arma::mat x(nr, nc);
  for (uint_t i = 0; i < nr; i++) {
    for (uint_t j = 0; j < nc; j++) {
      x(i, j) = std::pow(a, b(i, j));
    }
  }
  return x;
}
inline arma::mat flex_pow(const arma::mat& a, const double& b) {
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
inline arma::vec pnorm_cpp(const arma::vec& values, const bool& lower_tail) {
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
 Make matrices of coefficient estimates and standard errors, and matrix of covariances.
 */
inline void make_B_B_cov(arma::mat& B, arma::mat& B_cov, arma::vec& B0,
                         const arma::mat& iV,
                         const arma::mat& UU,
                         const arma::mat& X,
                         const std::vector<arma::mat>& U) {
  
  arma::mat mean_sd_X(X.n_cols, 2);
  mean_sd_X.col(0) = arma::conv_to<arma::vec>::from(arma::mean(X));
  mean_sd_X.col(1) = arma::conv_to<arma::vec>::from(arma::stddev(X));
  std::vector<arma::vec> sd_U(U.size());
  for (uint_t i = 0; i < U.size(); i++) {
    if (U[i].n_cols > 0) sd_U[i] = arma::conv_to<arma::vec>::from(arma::stddev(U[i]));
  }
  
  arma::vec sd_vec(UU.n_cols, arma::fill::zeros);
  
  for (uint_t counter = 0, i = 0; i < X.n_cols; counter++, i++) {
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



/*
 Returning useful error message if choleski decomposition fails:
 */
inline void safe_chol(arma::mat& L, std::string task) {
  try {
    L = arma::chol(L);
  } catch(const std::runtime_error& re) {
    std::string err_msg = static_cast<std::string>(re.what());
    if (err_msg == "chol(): decomposition failed") {
      std::string err_msg_out = "Choleski decomposition failed during " + task + ". ";
      err_msg_out += "Changing the `constrain_d` argument to `TRUE`, and ";
      err_msg_out += "using a different algorithm (`method` argument) can remedy this.";
      throw(Rcpp::exception(err_msg_out.c_str(), false));
    } else {
      std::string err_msg_out = "Runtime error: \n" + err_msg;
      throw(Rcpp::exception(err_msg_out.c_str(), false));
    }
  } catch(const std::exception& ex) {
    std::string err_msg = static_cast<std::string>(ex.what());
    stop("Error occurred: \n" + err_msg);
  } catch(...) {
    stop("Unknown failure occurred.");
  }
  return;
}



#endif
