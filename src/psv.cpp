// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
IntegerVector which2(const LogicalVector x){
  int n = x.size();
  int m = sum(x);
  if(m > n) stop("More true values than the length of the vector");
  IntegerVector xx(m);
  int cnt = 0;
  for(int i = 0; i < n; i++){
    if(x[i]){
      xx[cnt] = i;
      ++cnt;
    }
  }
  return xx;
}

// [[Rcpp::export]]
DataFrame psv_cpp(const NumericMatrix& comm, 
                  const arma::mat& Cmatrix,
                  const bool compute_var){
  int nlocations = comm.nrow();
  int nspecies = comm.ncol();
  NumericVector PSVs(nlocations); // to hold results
  NumericVector SR(nlocations); // to hold results
  for(int i = 0; i < nlocations; ++i){
    LogicalVector index = (comm.row(i) > 0);
    IntegerVector iindex = which2(index);
    int nsp = sum(index);
    double psv;
    if(nsp > 1){
      arma::uvec iindex_arma = as<arma::uvec>(iindex);
      arma::mat cm = Cmatrix.submat(iindex_arma, iindex_arma);
      psv = (nsp * trace(cm) - accu(cm)) / (nsp * (nsp - 1));
    } else {
      psv = NA_REAL;
    }
    PSVs[i] = psv;
    SR[i] = nsp;
  }
  
  NumericVector PSVary(nspecies - 1);
  NumericVector PSVvar(nspecies);
  if(compute_var && nspecies > 1){
    arma::mat cdig = mat(nspecies, nspecies);
    cdig.eye();
    arma::mat X = Cmatrix - (accu(Cmatrix - cdig))/(nspecies * (nspecies - 1));
    X.diag().zeros();
    double SS1 = accu(X % X) / 2;
    NumericMatrix ss2(nspecies);
    for(int i = 0; i < (nspecies - 1); i++){
      double sumi = sum(X.row(i));
      for(int j = i + 1; j < nspecies; j++){
        ss2(i, j) = X(i, j) * (sumi - X(i, j));
      }
    }
    double SS2 = sum(ss2);
    double SS3 = (-1 * SS1) - SS2;
    double S1 = SS1 * 2/(nspecies * (nspecies - 1));
    double S2 = SS2 * 2/(nspecies * (nspecies - 1) * (nspecies - 2));
    double S3;
    if (nspecies == 3) {
      S3 = 0;
    } else {
      S3 = SS3 * 2/(nspecies * (nspecies - 1) * (nspecies - 2) * (nspecies - 3));
    }
    
    for (int ni = 1; ni < (nspecies - 1); ni++) {
      double nii = ni + 1; // need to be double, not int
      PSVary[ni - 1] = 2/(nii * (nii - 1)) * (S1 + (nii - 2) * S2 + (nii - 2) * (nii - 3) * S3);
    }
    
    for(int g = 0; g < nlocations; g++){
      if(SR[g] > 1){
        PSVvar[g] = PSVary[(SR[g] - 2)];
      } else {
        PSVvar[g] = NA_REAL;
      }
    }
  }
  
  DataFrame PSVout = DataFrame::create( Named("PSVs") = PSVs , 
                                        Named("SR") = SR,
                                        Named("vars") = PSVvar);
  return PSVout;
}

/*** R
# psv_cpp(as.matrix(comm), Cmatrix, T)
# which(comm[2,] > 0)
# which2(comm[2,] > 0)
 */