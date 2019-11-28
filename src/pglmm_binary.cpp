// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List pglmm_iV_logdetV_cpp(NumericVector par, arma::vec mu,
                                const arma::sp_mat& Zt, const arma::sp_mat& St, 
                                const List& nested, bool logdet,
                                const std::string family, arma::vec totalSize){
  int q_nonNested = St.n_rows;
  arma::sp_mat Ut;
  arma::sp_mat U;
  if(q_nonNested > 0){
    IntegerVector idx = seq_len(q_nonNested) - 1; // c++ starts with 0
    // uvec idx_uvec = as<uvec>(idx);
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
    sn = as<NumericVector>(wrap(sn1)); // no need to declare type again
  } 
  
  arma::sp_mat iV0;
  arma::mat Ishort_Ut_iA_U;
  double logdetV;
  double signV;
  double logdetiA;
  double signiA;
  if (q_Nested == 0){ // then q_nonNested will not be 0, otherwise, no random terms
    arma::vec pq = totalSize % mu % (1 - mu);
    arma::sp_mat iA;
    if(family == "binomial") iA = sp_mat(diagmat(pq));
    if(family == "poisson") iA = sp_mat(diagmat(mu));
    arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
    arma::sp_mat Ut_iA_U = Ut * iA * U;
    // Woodbury identity
    Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
    arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
    iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
    arma::mat iV(iV0); // convert to dense matrix
    if(logdet){
      log_det(logdetV, signV, Ishort_Ut_iA_U); 
      log_det(logdetiA, signiA, mat(iA)); 
      logdetV = logdetV - logdetiA;
      NumericVector logdetV1 = NumericVector::create(logdetV);
      if(any(is_infinite(logdetV1))){
        arma::mat lgm = chol(Ishort_Ut_iA_U);
        logdetV = 2 * sum(log(lgm.diag())) - logdetiA;
      }
    }
  } else {
    arma::vec pq = 1 / (totalSize % mu % (1 - mu));
    arma::sp_mat A;
    if(family == "binomial") A = sp_mat(diagmat(pq));
    if(family == "poisson") A = sp_mat(diagmat(1 / mu));
    if (q_Nested == 1){
      double snj = pow(sn[0], 2);
      sp_mat nj = nested[0];
      A = A + snj * nj;
    } else {
      for (int j = 0; j < q_Nested; j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
    arma::mat A1(A);
    arma::sp_mat iA = sp_mat(inv(A1));
    // Rcout << iA << " " ;
    if(q_nonNested > 0){
      arma::sp_mat Ishort = sp_mat(Ut.n_rows, Ut.n_rows); Ishort.eye();
      arma::sp_mat Ut_iA_U = Ut * iA * U;
      Ishort_Ut_iA_U = mat(Ishort + Ut_iA_U);
      arma::mat i_Ishort_Ut_iA_U = inv(Ishort_Ut_iA_U);
      iV0 = iA - iA * U * sp_mat(i_Ishort_Ut_iA_U) * Ut * iA;
    } else {
      iV0 = iA;
    }
    
    arma::mat iV(iV0); // convert to dense matrix
    if(logdet){
      log_det(logdetV, signV, iV); 
      logdetV = -1 * logdetV;
      NumericVector logdetV1 = NumericVector::create(logdetV);
      if(any(is_infinite(logdetV1))){
        arma::mat lgm = chol(iV);
        logdetV = -2 * sum(log(lgm.diag()));
      }
    }
  }
  
  if(logdet){
    return List::create(
      _["iV"] = iV0,
      _["logdetV"] = logdetV
    );
  } else {
    return List::create(_["iV"] = iV0);
  }
}

// [[Rcpp::export]]
arma::sp_mat pglmm_V(NumericVector par, const arma::sp_mat& Zt, 
                           const arma::sp_mat& St, arma::vec mu, 
                           const List& nested, bool missing_mu,
                           const std::string family, arma::vec totalSize){
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
    sn = as<NumericVector>(wrap(sn1)); // no need to declare type again
  } 
  
  arma::mat iW;
  if(missing_mu){
    iW = mat(Zt.n_cols, Zt.n_cols, fill::zeros);
  } else {
    arma::vec pq = 1 / (totalSize % mu % (1 - mu));
    if(family == "binomial") iW = diagmat(pq);
    if(family == "poisson") iW = diagmat(1 / mu);
  }
  
  arma::sp_mat A = sp_mat(iW);
  if(q_Nested > 0){
    if (q_Nested == 1){
      double snj = pow(sn[0], 2);
      sp_mat nj = nested[0];
      A = A + snj * nj;
    } else {
      for (int j = 0; j < q_Nested; j++) {
        double snj = pow(sn[j], 2);
        sp_mat nj = nested[j];
        A = A + snj * nj;
      }
    }
  }
  
  arma::sp_mat V;
  if(q_nonNested > 0){
    V = A + U * Ut;
  } else {
    V = A;
  }
  
  return V;
}

// [[Rcpp::export]]
double pglmm_LL_cpp(NumericVector par, const arma::vec& H,
                          const arma::mat& X, const arma::sp_mat& Zt, 
                          const arma::sp_mat& St, const arma::vec& mu, 
                          const List& nested, bool REML, bool verbose,
                          const std::string family, arma::vec totalSize){
  par = abs(par);
  // unsigned int n = H.n_rows;
  List iVdet = pglmm_iV_logdetV_cpp(par, mu, Zt, St, nested, true, family, totalSize);
  sp_mat iV0 = as<sp_mat>(iVdet["iV"]);
  mat iV = mat(iV0);
  double logdetV = iVdet["logdetV"];
  double LL;
  if (REML) {
    double logdetL, signL;
    log_det(logdetL, signL, trans(X) * iV * X);
    LL = 0.5 * (logdetV + as_scalar(trans(H) * iV * H) + logdetL);
  } else {
    LL = 0.5 * (logdetV + as_scalar(trans(H) * iV * H));
  }
  
  if (verbose) Rcout << LL << " " << par << std::endl;
  
  return LL;
}


// [[Rcpp::export]]
List pglmm_internal_cpp(const arma::mat& X, const arma::vec& Y,
                               const arma::sp_mat& Zt, const arma::sp_mat& St,
                               const List& nested, const bool REML, const bool verbose,
                               const int n, const int p, const int q, const int maxit, 
                               const double reltol, const double tol_pql, const double maxit_pql,
                               const std::string optimizer, arma::mat B_init, arma::vec ss,
                               const std::string family, arma::vec totalSize){
  Rcpp::checkUserInterrupt();
  mat B = B_init;
  mat b(n, 1, fill::zeros);
  mat beta = join_vert(B, b);
  vec mu;
  if(family == "binomial") mu = arma::exp(X * B) / (1 + arma::exp(X * B));
  if(family == "poisson") mu = arma::exp(X * B);
  mat ix(n, n, fill::eye);
  mat XX = join_horiz(X, ix);
  
  vec est_ss = ss;
  vec est_B = B;
  vec oldest_ss(size(ss));
  oldest_ss.fill(1000000.0);
  mat oldest_B(size(B));
  oldest_B.fill(1000000.0);
  
  unsigned int iteration = 0, iteration_m;
  double tol_pql2 = pow(tol_pql, 2);
  double LL;
  
  NumericVector ss0 = wrap(ss); // to work with other functions
  vec Z, H, niter;
  int convcode;
  mat iV;
  
  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];
  Rcpp::Environment nloptr_pkg = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptr_pkg["nloptr"];
  Rcpp::Environment phyr_pkg = Rcpp::Environment::namespace_env("phyr");
  Rcpp::Function pglmm_LL_cpp2 = phyr_pkg["pglmm_LL_cpp"];
  
  while((as_scalar(trans(est_ss - oldest_ss) * (est_ss - oldest_ss)) > tol_pql2 ||
        as_scalar(trans(est_B - oldest_B) * (est_B - oldest_B)) > tol_pql2) &&
        iteration <= maxit_pql){
    oldest_ss = est_ss;
    oldest_B = est_B;
    vec est_B_m = B;
    mat oldest_B_m(size(est_B));
    oldest_B_m.fill(1000000.0);
    iteration_m = 0;
    Rcpp::checkUserInterrupt();
    
    // mean component
    while(as_scalar(trans(est_B_m - oldest_B_m) * (est_B_m - oldest_B_m)) > tol_pql2 &&
          iteration_m <= maxit_pql){
      // Rcpp::checkUserInterrupt();
      oldest_B_m = est_B_m;
      List iv = pglmm_iV_logdetV_cpp(ss0, mu, Zt, St, nested, false, family, totalSize);
      sp_mat iV0 = iv["iV"];
      if(family == "binomial") Z = X * B + b + (Y/totalSize - mu)/(mu % (1 - mu));
      if(family == "poisson") Z = X * B + b + (Y - mu)/mu;
      
      iV = mat(iV0); // convert to dense matrix
      arma::mat denom = trans(X) * iV * X;
      arma::mat num = trans(X) * iV * Z;
      B = solve(denom, num);
      
      sp_mat V = pglmm_V(ss0, Zt, St, mu, nested, false, family, totalSize);
      vec diav = vectorise(1/(totalSize % mu % (1 - mu)));
      sp_mat iW;
      if(family == "binomial") iW = sp_mat(diagmat(diav));
      if(family == "poisson") iW = sp_mat(diagmat(1/mu));
      sp_mat C = V - iW;
      b = mat(C) * mat(iV) * (Z - X * B);
      beta = join_vert(B, b);
      if(family == "binomial") mu = arma::exp(XX * beta) / (1 + arma::exp(XX * beta));
      if(family == "poisson") mu = arma::exp(XX * beta);
      
      est_B_m = B;
      if(verbose) Rcout << "mean part: " << iteration_m << " " << trans(B) << std::endl;
      ++iteration_m;
      // Rcout << "mean part: " << iteration_m << " " << trans(B) << std::endl;
      // Rcout << "            denom: " << denom << endl;
      // Rcout << "            num: " << num << endl;
      if(B.has_nan()) Rcpp::stop("Estimation of B failed. Check for lack of variation in Y. You could try with a smaller s2.init, but this might not help.");
    } // end while for mean
    
    // variance component
    if(family == "binomial") Z = X * B + b + (Y/totalSize - mu)/(mu % (1 - mu)); // B, b, mu all updated
    if(family == "poisson") Z = X * B + b + (Y - mu)/mu;
    H = Z - X * B;
   
    Rcpp::List opt;
    if(optimizer == "Nelder-Mead"){
      if(q > 1){
        opt = optim(_["par"] = ss0,
                    _["fn"] = pglmm_LL_cpp2,
                    _["H"] = H, _["X"] = X, _["Zt"] = Zt,
                    _["St"] = St, _["mu"] = mu, _["nested"] = nested,
                      _["REML"] = REML, _["verbose"] = verbose,
                      _["family"] = family, _["totalSize"] = totalSize,
                      _["method"] = "Nelder-Mead",
                      _["control"] = List::create(_["maxit"] = maxit, _["reltol"] = reltol));
      } else {
        // opt = optim(_["par"] = ss0,
        //             _["fn"] = Rcpp::InternalFunction(&plmm_binary_LL_cpp),
        //             _["X"] = X, _["H"] = H, _["Zt"] = Zt,
        //             _["St"] = St, _["nested"] = nested,
        //             _["mu"] = mu, _["REML"] = REML, _["verbose"] = verbose,
        //               _["method"] = "L-BFGS-B",
        //               _["control"] = List::create(_["maxit"] = maxit, _["factr"] = reltol));
        Rcpp::stop("With only 1 random term and cpp = TRUE, phyr cannot run the optimization yet. \n \
                     Set optimizer to other options, e.g. nelder-mead-nlopt and re-run it. \n \
                     Or you can turn cpp off with cpp = FALSE and re-run it.");
      }
    } else {
      std::string nlopt_algor;
      if (optimizer == "bobyqa") nlopt_algor = "NLOPT_LN_BOBYQA";
      if (optimizer == "nelder-mead-nlopt") nlopt_algor = "NLOPT_LN_NELDERMEAD";
      if (optimizer == "subplex") nlopt_algor = "NLOPT_LN_SBPLX";
      List opts = List::create(_["algorithm"] = nlopt_algor,
                               _["ftol_rel"] = reltol, _["ftol_abs"] = reltol,
                               _["xtol_rel"] = 0.0001,
                               _["maxeval"] = maxit);
      List S0 = nloptr(_["x0"] = ss0,
                       _["eval_f"] = pglmm_LL_cpp2,
                       _["opts"] = opts, _["H"] = H, _["X"] = X, _["Zt"] = Zt,
                       _["St"] = St, _["mu"] = mu, _["nested"] = nested,
                         _["REML"] = REML, _["verbose"] = verbose,
                         _["family"] = family, _["totalSize"] = totalSize);
      opt = List::create(_["par"] = S0["solution"], _["value"] = S0["objective"], 
                         _["counts"] = S0["iterations"], _["convergence"] = S0["status"], 
                         _["message"] = S0["message"]);
    }
      
    arma::vec par_opt0 = abs(as<arma::vec>(opt["par"]));
    ss0 = wrap(par_opt0);
    LL = as_scalar(as<double>(opt["value"]));
    convcode = as<int>(opt["convergence"]);
    niter = as<arma::vec>(opt["counts"]);
    
    est_ss = par_opt0;
    est_B = B;
    ++iteration;
    if(verbose) Rcout << "var part: " << iteration << " " << LL << std::endl;
    // Rcout << "var part: " << iteration << " " << LL << " " << ss0 << std::endl;
    // } // end opt
  } // end while
  
  List out = List::create(
    _["B"] = B, _["ss"] = ss0, 
    _["iV"] = iV, _["mu"] = mu, _["H"] = H,
      _["convcode"] = convcode,
      _["niter"] = niter,
      _["LL"] = LL
  );
  
  return out;
}

// [[Rcpp::export]]
int sexp_type(SEXP x){ 
  return TYPEOF(x); 
}


/*** R
# library(Matrix)
# St = as(matrix(NA, 0, 0), "dgTMatrix")
# Zt = as(matrix(NA, 0, 0), "dgTMatrix")
# sexp_type(nested)
# plmm_binary_iV_logdetV_cpp(ss, mu, Zt, St, nested, F)
# plmm_binary_V(ss, Zt, St, mu, nested, F)
# plmm_binary_LL_cpp(ss, mu, X, Zt, St, mu, nested, T, T)
# internal_res = pglmm_binary_internal_cpp(X = X, Y = Y, Zt = Zt, St = St, 
#                                          nested = nested, REML = REML, verbose = verbose, 
#                                          n = n, p = p, q = q, maxit = maxit, 
#                                          reltol = reltol, tol_pql = tol.pql, 
#                                          maxit_pql = maxit.pql, optimizer = "bobyqa", 
#                                          B_init = B.init, ss = ss)
# phyr:::pglmm_binary_internal_cpp(X = X, Y = Y, Zt = Zt, St = St, 
#                           nested = nested, REML = REML, verbose = verbose, 
#                           n = n, p = p, q = q, maxit = maxit, 
#                           reltol = reltol, tol_pql = tol.pql, 
#                           maxit_pql = maxit.pql, optimizer = optimizer, 
#                           B_init = B.init, ss = ss)
# opt <- optim(fn = plmm.binary.LL, par = ss, H = H, X = X, Zt = Zt, St = St,
#              mu = mu, nested = nested, REML = REML, verbose = verbose, 
#              method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
# opt2 <- optim(fn = plmm_binary_LL_cpp, par = ss, H = as.matrix(H), X = X, Zt = Zt, St = St,
#              mu = mu, nested = nested, REML = REML, verbose = T,
#              method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))

# plmm_binary_LL_cpp(ss, as.matrix(H), X, Zt,  St, mu,  nested, REML, verbose)
# plmm_binary_LL_cpp(c(0.5, 0.5, 0.5, 0.5), as.matrix(H), X, Zt, St, mu, nested, REML = F, T)
*/