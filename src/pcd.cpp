// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
NumericVector predict_cpp(int n, const arma::vec& nsr, int reps, const arma::mat& V){
  int n_unique = nsr.size();
  NumericVector SSii(n_unique); // to save results
  int n1 = 2; // the number of n1 does not matter
  IntegerVector n1_idx = seq_len(n1) - 1;

  for(int n2 = 0; n2 < n_unique; n2++){
    NumericVector temp(reps);
    IntegerVector n2_idx = seq_len(nsr[n2]) - 1;
    for(int t = 0; t < reps; t++){
      // set_seed(t + 1);
      IntegerVector rp = sample(n, n);
      IntegerVector pick1 = rp[n1_idx];
      IntegerVector pick1_1 = pick1 - 1;
      arma::uvec pick1_arma = as<arma::uvec>(pick1_1); // convert to armadillo
      
      // set_seed(t + 10000);
      IntegerVector rp2 = sample(n, n);
      IntegerVector pick2 = rp2[n2_idx];
      IntegerVector pick2_1 = pick2 - 1;
      arma::uvec pick2_arma = as<arma::uvec>(pick2_1);
      
      arma::mat C11 = V.submat(pick1_arma, pick1_arma);
      arma::mat C22 = V.submat(pick2_arma, pick2_arma);
      arma::mat C12 = V.submat(pick1_arma, pick2_arma);
      
      arma::mat invC22 = inv(C22);
      arma::mat S11 = C11 - C12 * invC22 * trans(C12);
      double SS11 = (n1 * trace(S11) - accu(S11)) / (n1 * (n1 - 1));
      temp[t] = SS11;
    }
    SSii[n2] = mean(temp);
  }
  return SSii;
}

// [[Rcpp::export]]
List pcd2_loop(arma::vec SSii, arma::vec nsr, double SCii, const arma::mat& comm, 
               const arma::mat& V, int nsp_pool, bool verbose){
  int m = comm.n_rows;
  int n = comm.n_cols;
  NumericMatrix PCD(m, m);
  NumericMatrix PCDc(m, m);
  NumericMatrix PCDp(m, m);
  NumericMatrix D_pairwise(m, m);
  NumericMatrix dsor_pairwise(m, m);
  IntegerVector nn = seq_len(n);
  vec nna = as<vec>(nn);
  
  for(int i = 0; i < (m - 1); i++){
    if (verbose) {Rcout << i + 1 << " " ;}
    for(int j = i + 1; j < m; j++){
      // Rcout << "Rows: " << j + 1 << std::endl ;
      uvec c_i = find(comm.row(i) == 1);
      uvec c_j = find(comm.row(j) == 1);
      vec pick1 = nna.elem(c_i);
      vec pick2 = nna.elem(c_j);
      IntegerVector pick1_cpp = wrap(pick1);
      IntegerVector pick2_cpp = wrap(pick2);
      IntegerVector pick_inter = intersect(pick1_cpp, pick2_cpp);
      
      vec pick12 = join_cols(pick1, pick2); // how to convert to a uvec?
      IntegerVector pick12_int = wrap(pick12);
      IntegerVector pick12_intc = pick12_int - 1; // c++ starts at 0
      uvec pick12_uvec = as<uvec>(pick12_intc);
      
      int n1 = pick1.size();
      int n2 = pick2.size();
      
      mat C = V.submat(pick12_uvec, pick12_uvec);
      mat C11 = C.submat(0, 0, n1 - 1, n1 - 1);
      mat C22 = C.submat(n1, n1, n1+n2-1, n1+n2-1);
      mat C12 = C.submat(0, n1, n1-1, n1+n2-1);
      
      mat invC11 = inv(C11);
      mat S22 = C22 - trans(C12) * invC11 * C12; // n2 by n2 matrix
      
      mat invC22 = inv(C22);
      mat S11 = C11 - C12 * invC22 * trans(C12); // n2 by n2 matrix
      
      double SC11;
      double SS11;
      double SC22;
      double SS22;
      
      if(n1 > 1){
        SC11 = (n1 * trace(C11) - accu(C11)) / double(n1 * (n1 - 1));
        SS11 = (n1 * trace(S11) - accu(S11)) / double(n1 * (n1 - 1));
      } else {
        SC11 = (n1 * trace(C11) - accu(C11)) / double(n1 * n1);
        SS11 = (n1 * trace(S11) - accu(S11)) / double(n1 * n1);
      }
      
      if(n2 > 1){
        SC22 = (n2 * trace(C22) - accu(C22)) / double(n2 * (n2 - 1));
        SS22 = (n2 * trace(S22) - accu(S22)) / double(n2 * (n2 - 1));
      } else {
        SC22 = (n2 * trace(C22) - accu(C22)) / double(n2 * n2);
        SS22 = (n2 * trace(S22) - accu(S22)) / double(n2 * n2);
      }
      double D = (n1 * SS11 + n2 * SS22) / double(n1 * SC11 + n2 * SC22);
      double dsor = 1 - 2 * pick_inter.size() / double(n1 + n2);
      uvec which_n2 = find(nsr == n2) ;
      uvec which_n1 = find(nsr == n1) ;
      double ssii_n2 = as_scalar(SSii.elem(which_n2));
      double ssii_n1 = as_scalar(SSii.elem(which_n1));
      double pred_D = (n1 * ssii_n2 + n2 * ssii_n1) / (n1 * SCii + n2 * SCii);
      double pred_dsor = 1 - 2 * n1 * n2 / double((n1 + n2) * nsp_pool);
      
      double PCDij = D/pred_D;
      double PCDcij = dsor/pred_dsor;
      double PCDpij = PCDij/PCDcij;
      PCD(i, j) = PCDij;
      PCDc(i, j) = PCDcij;
      PCDp(i, j) = PCDpij;
      D_pairwise(i, j) = D;
      dsor_pairwise(i, j) = dsor;
    }
  }
  
  return List::create(
    _["PCD"] = PCD,
    _["PCDc"] = PCDc,
    _["PCDp"] = PCDp,
    _["D_pairwise"] = D_pairwise,
    _["dsor_pairwise"] = dsor_pairwise
  );
}