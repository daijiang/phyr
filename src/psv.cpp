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
arma::mat vcv_loop(NumericVector& xx, const int& n, const IntegerVector& e1, 
                   const IntegerVector& e2, const NumericVector& EL, 
                   const List& pp, bool corr //, const CharacterVector& sp
  ){
  arma::mat vcv(n, n, fill::zeros);
  int len_e1 = e1.size();
  for(int i = (len_e1 - 1); i > (-1); i--){
    double var_cur_node = xx[e1[i] - 1];
    xx[e2[i] - 1] = var_cur_node + EL[i];
    int j = (i - 1);
    while((e1[j] == e1[i]) && (j > (-1))){
      IntegerVector left, right;
      if(e2[j] > n){
        left = pp[e2[j] - n - 1];
      } else {
        left = e2[j];
      }
      if(e2[i] > n){
        right = pp[e2[i] - n - 1];
      } else {
        right = e2[i];
      }
      left = left -1 ;
      right = right - 1;
      arma::uvec left2, right2;
      left2 = as<arma::uvec>(left);
      right2 = as<arma::uvec>(right);
      vcv.submat(left2, right2).fill(var_cur_node);
      vcv.submat(right2, left2).fill(var_cur_node);
      j = j - 1;
    }
  }
  arma::vec xx2 = as<arma::vec>(xx);
  vcv.diag() = xx2.subvec(0, n - 1);
  if(corr){
    arma::vec Is = sqrt(1 / vcv.diag());
    vcv.each_col() %= Is;
    vcv.each_row() %= trans(Is);
    vcv.diag().ones();
  }
  // NumericMatrix vcvv = wrap(vcv);
  // rownames(vcvv) = sp;
  // colnames(vcvv) = sp;
  return vcv;
}

// [[Rcpp::export]]
void cov2cor_cpp(arma::mat& vcv){
  arma::vec Is = sqrt(1 / vcv.diag());
  vcv.each_col() %= Is;
  vcv.each_row() %= trans(Is);
  vcv.diag().ones();
  return ;
}

// [[Rcpp::export]]
NumericVector pse_cpp(const NumericMatrix& comm, const arma::mat& Cmatrix){
  int nlocations = comm.nrow();
  // int nspecies = comm.ncol();
  NumericVector PSEs(nlocations); // to hold results
  // NumericVector SR(nlocations); // to hold results
  for(int i = 0; i < nlocations; ++i){
    Rcpp::checkUserInterrupt();
    LogicalVector index = (comm.row(i) > 0);
    IntegerVector iindex = which2(index);
    int nsp = sum(index);
    double pse;
    if(nsp > 1){
      arma::uvec iindex_arma = as<arma::uvec>(iindex);
      arma::mat C = Cmatrix.submat(iindex_arma, iindex_arma);
      NumericVector M0 = comm(i, _);
      double N = sum(M0);
      NumericVector M1 = M0[iindex];
      arma::colvec M = as<arma::colvec>(M1);
      double mbar = mean(M1);
      double dem = pow(N, 2) - N * mbar;
      double num = as_scalar(N * trans(C.diag()) * M - trans(M) * C * M);
      pse = (num / dem);
      // Rcout << pse << " ";
    } else {
      pse = NA_REAL;
    }
    PSEs[i] = pse;
  }
  return PSEs;
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
  
  Rcpp::checkUserInterrupt();
  
  NumericVector PSVvar(nlocations);
  if(compute_var && nspecies > 1){
    arma::mat cdig = mat(nspecies, nspecies);
    cdig.eye();
    arma::mat X = Cmatrix - (accu(Cmatrix - cdig))/(nspecies * (nspecies - 1));
    X.diag().zeros();
    double SS1 = accu(X % X) / 2;
    // Rcout << "SS1 = " << SS1 << endl;
    NumericMatrix ss2(nspecies);
    for(int i = 0; i < (nspecies - 1); i++){
      double sumi = sum(X.row(i));
      for(int j = i + 1; j < nspecies; j++){
        ss2(i, j) = X(i, j) * (sumi - X(i, j));
      }
    }
    double SS2 = sum(ss2);
    double SS3 = (-1 * SS1) - SS2;
    // Rcout << "SS2 = " << SS2 << endl;
    // Rcout << "SS3 = " << SS3 << endl;
    // Rcout << "nspecies * (nspecies - 1) = " << nspecies * (nspecies - 1) << endl;
    // Rcout << "(nspecies * (nspecies - 1) * (nspecies - 2)) = " << double(nspecies * (nspecies - 1) * (nspecies - 2)) << endl;
    double S1 = SS1 * 2.0/(nspecies * (nspecies - 1.0));
    double S2 = SS2 * 2.0/(nspecies * (nspecies - 1.0) * (nspecies - 2.0)); // 1.0 not 1 !!
    double S3;
    if (nspecies == 3) {
      S3 = 0;
    } else {
      S3 = SS3 * 2.0/(nspecies * (nspecies - 1.0) * (nspecies - 2.0) * (nspecies - 3.0));
    }
    // Rcout << "S1 = " << S1 << endl;
    // Rcout << "S2 = " << S2 << endl;
    // Rcout << "S3 = " << S3 << endl;
    // Rcout << "nsp = " << nspecies << endl;
    Rcpp::checkUserInterrupt();
    NumericVector PSVary(nspecies - 1);
    for (int ni = 1; ni < (nspecies - 1); ni++) {
      double nii = ni + 1; // need to be double, not int
      PSVary[ni - 1] = 2.0/(nii * (nii - 1)) * (S1 + (nii - 2) * S2 + (nii - 2) * (nii - 3) * S3);
    }
    
    for(int g = 0; g < nlocations; g++){
      if(SR[g] > 1){
        PSVvar[g] = PSVary[(SR[g] - 2)];
      } else {
        PSVvar[g] = NA_REAL;
      }
    }
  }

  return DataFrame::create(
    _["PSVs"] = PSVs, 
    _["SR"] = SR,
    _["vars"] = PSVvar
  );
}

/*** R
# nspp = 20
# nsite = 30
# comm_sim = matrix(rbinom(nspp * nsite, size = 1, prob = 0.6), nrow = nsite, ncol = nspp)
# row.names(comm_sim) = paste0("site_", 1:nsite)
# colnames(comm_sim) = paste0("t", 1:nspp)
# tree_sim = ape::rtree(n = nspp)
# comm_sim = comm_sim[, tree_sim$tip.label]
# cmx = ape::vcv(tree_sim)
# test = psv_cpp(comm_sim, cmx, TRUE)
# test
# psv_cpp(comm_sim, cmx, F)
# tst = psv_cpp(comm, Cmatrix, T)
# psv_cpp(comm, Cmatrix, T)
*/