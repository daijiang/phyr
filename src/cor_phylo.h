// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef __PHYR_COR_PHYLO_H
#define __PHYR_COR_PHYLO_H

#include <RcppArmadillo.h>
#include <numeric>
#include <cmath>
#include <vector>
#include <math.h>


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
class LogLikInfo {
public:
  arma::vec par0;  // par to start with
  arma::mat XX;
  arma::mat UU;
  arma::mat MM;
  arma::mat Vphy;
  arma::mat tau;
  bool REML;
  bool no_corr;
  bool constrain_d;
  double lower_d;
  bool verbose;
  double rcond_threshold;
  uint_t iters;
  arma::vec min_par; // par for minimum LL
  double LL;
  int convcode;
  
  LogLikInfo() {}
  LogLikInfo(const arma::mat& X,
          const std::vector<arma::mat>& U,
          const arma::mat& M,
          const arma::mat& Vphy_,
          const bool& REML_,
          const bool& no_corr_,
          const bool& constrain_d_,
          const double& lower_d_,
          const bool& verbose_,
          const double& rcond_threshold_);
  // Used in bootstrapping
  LogLikInfo(const arma::mat& X,
          const std::vector<arma::mat>& U,
          const arma::mat& M,
          XPtr<LogLikInfo> other);
  // Copy constructor
  LogLikInfo(const LogLikInfo& ll_info2) {
    par0 = ll_info2.par0;
    XX = ll_info2.XX;
    UU = ll_info2.UU;
    MM = ll_info2.MM;
    Vphy = ll_info2.Vphy;
    tau = ll_info2.tau;
    REML = ll_info2.REML;
    no_corr = ll_info2.no_corr;
    constrain_d = ll_info2.constrain_d;
    lower_d = ll_info2.lower_d;
    verbose = ll_info2.verbose;
    rcond_threshold = ll_info2.rcond_threshold;
    iters = ll_info2.iters;
    min_par = ll_info2.min_par;
    LL = ll_info2.LL;
    convcode = ll_info2.convcode;
  }
  
};



// Results from bootstrapping

class BootResults {
public:
  arma::cube corrs;
  arma::mat B0;
  arma::cube B_cov;
  arma::mat d;
  std::vector<arma::mat> out_mats;
  std::vector<uint_t> out_inds;
  std::vector<int> out_codes;

  BootResults(const uint_t& p, const uint_t& B_rows, const uint_t& n_reps) 
    : corrs(p, p, n_reps, arma::fill::zeros), 
      B0(B_rows, n_reps, arma::fill::zeros), 
      B_cov(B_rows, B_rows, n_reps, arma::fill::zeros),
      d(p, n_reps, arma::fill::zeros), 
      out_mats(), out_inds(), out_codes() {};

  // Insert values into a BootResults object
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
class BootMats {
public:
  // original input matrices
  const arma::mat X;
  const std::vector<arma::mat> U;
  const arma::mat M;
  arma::mat X_new;
  
  BootMats(const arma::mat& X_, const std::vector<arma::mat>& U_,
            const arma::mat& M_,
            const arma::mat& B_, const arma::vec& d_, XPtr<LogLikInfo> ll_info);
  
  XPtr<LogLikInfo> iterate(XPtr<LogLikInfo> ll_info);
  
  void one_boot(XPtr<LogLikInfo> ll_info, BootResults& br,
                const uint_t& i, const double& rel_tol, const int& max_iter,
                const std::string& method, const std::string& keep_boots,
                const std::vector<double>& sann);
  
  
private:
  arma::mat iD;
  arma::mat X_pred;

  // Method for returning bootstrapped data
  void boot_data(XPtr<LogLikInfo> ll_info, BootResults& br, const uint_t& i);

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


// pnorm for standard normal (i.e., ~ N(0,1))
inline arma::vec pnorm_cpp(const arma::vec& values, const bool& lower_tail) {
  
  arma::vec out = -1 * values * M_SQRT1_2;
  for (uint_t ii = 0; ii < out.n_elem; ii++) out(ii) = erfc(out(ii));
  out *= 0.5;
  
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


inline arma::vec make_par(const uint_t& p, const arma::mat& L, const bool& no_corr) {
  
  
  
  if (!no_corr) {
    
    uint_t par_size = (static_cast<double>(p) / 2.0) * (1 + p) + p;
    arma::vec par0(par_size);
    par0.fill(0.5);
    
    for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
      par0(arma::span(j, k)) = L(arma::span(i, p-1), i);
      j = k + 1;
      k += (p - i - 1);
    }
    
    return par0;
    
  }
  
  
  arma::vec par0(p * 2);
  par0.fill(0.5);
  arma::vec Ldiag = L.diag();
  for (uint_t i = 0; i < p; i++) par0(i) = Ldiag(i);
  
  return par0;
  
}


inline arma::mat make_L(const arma::vec& par, const uint_t& p) {
  
  arma::mat L(p, p, arma::fill::zeros);
  
  if (par.n_elem == static_cast<uint_t>((static_cast<double>(p) / 2) * (1 + p) + p)) {
    
    for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
      L(arma::span(i, p-1), i) = par(arma::span(j, k));
      j = k + 1;
      k += (p - i - 1);
    }
    
  } else if (par.n_elem == (2 * p)) {
    
    for (uint_t i = 0; i < p; i++) {
      L(i, i) = par(i);
    }
    
  } else {
    
    stop("\nINTERNAL ERROR: inappropriate length of `par` inside `make_L`");
    
  }
  
  return L;
  
}
inline arma::mat make_L(NumericVector par, const uint_t& p) {
  
  arma::mat L(p, p, arma::fill::zeros);
  
  if (par.size() == static_cast<int>((static_cast<double>(p) / 2) * (1 + p) + p)) {
    
    for (uint_t i = 0, j = 0, k = p - 1; i < p; i++) {
      for (uint_t l = 0; l < (k-j+1); l++) L(i+l, i) = par[j+l];
      j = k + 1;
      k += (p - i - 1);
    }
    
  } else if (par.size() == static_cast<int>(2 * p)) {
    
    for (uint_t i = 0; i < p; i++) {
      L(i, i) = par[i];
    }
    
  } else {
    
    stop("\nINTERNAL ERROR: inappropriate length of `par` inside `make_L`");
    
  }
  
  return L;
  
}




inline arma::vec make_d(NumericVector par, 
                        const uint_t& p,
                        const bool& constrain_d, 
                        const double& lower_d,
                        bool& return_max) {
  arma::vec d(p, arma::fill::zeros);
  return_max = false;
  uint_t size_ = par.size();
  if (constrain_d) {
    arma::vec logit_d(p, arma::fill::zeros);
    for (uint_t i = 0, j = (size_ - p); j < size_; i++, j++) {
      logit_d(i) = par[j];
    }
    /*  --------------------------------  */
    // In function `cor_phylo_LL`, `return_max = true` indicates to return a huge value
    if (arma::max(arma::abs(logit_d)) > 10) {
      return_max = true;
      return d;
    }
    /*  --------------------------------  */
    
    for (uint_t i = 0; i < p; i++) d(i) = 1/(1 + std::exp(-1 * logit_d(i)));
    // If you ever want to allow this to be changed:
    double upper_d = 1.0;
    d *= (upper_d - lower_d);
    d += lower_d;
  } else {
    d.set_size(p);
    for (uint_t i = 0, j = (size_ - p); j < size_; i++, j++) {
      d(i) = par[j];
    }
    d += lower_d;
    /*  --------------------------------  */
    if (arma::max(d) > 10) return_max = true;
    /*  --------------------------------  */
  }
  return d;
}
inline arma::vec make_d(const arma::vec& par, const uint_t& p,
                        const bool& constrain_d, const double& lower_d) {
  arma::vec d(p, arma::fill::zeros);
  if (constrain_d) {
    arma::vec logit_d = par.tail(p);
    for (uint_t i = 0; i < p; i++) d(i) = 1 / (1 + std::exp(-1 * logit_d(i)));
    // If you ever want to allow this to be changed:
    double upper_d = 1.0;
    d *= (upper_d - lower_d);
    d += lower_d;
  } else {
    d = par.tail(p);
    d += lower_d;
  }
  return d;
}


// OU transform
inline arma::mat make_C(const uint_t& n, const uint_t& p,
                        const arma::mat& tau, const arma::vec& d, 
                        const arma::mat& Vphy, const arma::mat& R) {
  
  
  arma::mat C(p * n, p * n, arma::fill::zeros);
  arma::mat Cd, w, x, y, z;
  for (uint_t i = 0; i < p; i++) {
    for (uint_t j = 0; j < p; j++) {
      x = flex_pow(d(i), tau);
      y = flex_pow(d(j), tau.t());
      w = flex_pow(d(i) * d(j), Vphy);
      z = 1 - w;
      Cd = x % y % z;
      Cd /= (1 - d(i) * d(j));
      Cd *= R(i,j);
      for (uint_t ii = n * i, kk = 0; ii < (i + 1) * n; ii++, kk++) {
        for (uint_t jj = n * j, ll = 0; jj < (j + 1) * n; jj++, ll++) {
          C(ii,jj) = Cd(kk,ll);
        }
      }
    }
  }
  
  return C;
}

inline arma::mat make_V(const arma::mat& C, const arma::mat& MM) {
  arma::vec MMvec = arma::vectorise(MM);
  arma::mat MMdiagmat = arma::diagmat(MMvec);
  arma::mat V(C.n_rows, C.n_cols, arma::fill::zeros);
  for (uint_t i = 0; i < C.n_rows; i++) {
    for (uint_t j = 0; j < C.n_cols; j++) {
      V(i,j) = C(i,j) + MMdiagmat(i,j);
    }
  }
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
  
  B_cov = arma::inv(UU.t() * iV * UU);
  B_cov = arma::diagmat(sd_vec) * B_cov * arma::diagmat(sd_vec);
  
  B.set_size(B0.n_elem, 4);
  B.col(0) = B0;  // Estimates
  B.col(1) = arma::diagvec(B_cov);
  for (uint_t i = 0; i < B.n_rows; i++) B(i,1) = std::sqrt(B(i,1));
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
