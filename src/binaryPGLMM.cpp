// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// 
// // we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"
// 
// // via the depends attribute we tell Rcpp to create hooks for
// // RcppArmadillo so that the build process will know what to do
// //
// // [[Rcpp::depends(RcppArmadillo)]]
// 
// using namespace Rcpp;
// using namespace arma;
// 
// // [[Rcpp::export]]
// double pglmm_reml_cpp(arma::vec par, const arma::mat& tinvW,
//                       const arma::mat& tH, const arma::mat& tVphy,
//                       const arma::mat& tX){
//   vec ss2 = abs(real(par));
//   arma::mat Cd = as_scalar(ss2) * tVphy;
//   double LL = pow(static_cast<double>(10), static_cast<double>(10));
//   arma::mat V = tinvW + Cd;
//   if(!V.has_inf()){
//     cx_vec eigval = eig_gen(V);
//     vec eigreal = real(eigval);
//     if(all(eigreal > 0.0)){
//       arma::mat invV = inv(V);
//       double logdetV;
//       double signV;
//       log_det(logdetV, signV, V);
//       NumericVector logdetV1 = NumericVector::create(logdetV);
//       if(any(is_infinite(logdetV1))){
//         arma::mat lgm = chol(V);
//         logdetV = 2 * sum(log(lgm.diag()));
//       }
//       arma::mat x = trans(tX) * invV * tX;
//       double logdetx;
//       double signx;
//       log_det(logdetx, signx, x);
//       LL = logdetV + as_scalar(trans(tH) * invV * tH) + logdetx;
//     }
//   }
//   return LL;
// }
// 
// // [[Rcpp::export]]
// List binpglmm_inter_while_cpp(
//     arma::mat est_B_m, 
//     arma::mat oldest_B_m, 
//     arma::mat B, 
//     const double& tol_pql, 
//     int iteration_m, 
//     const int& maxit_pql,
//     arma::vec mu, 
//     const arma::mat& C, 
//     int rcondflag, 
//     const arma::mat& B_init, 
//     const arma::mat& X, 
//     const arma::mat& XX,
//     const arma::mat& est_B, 
//     const arma::vec& y, 
//     int n, arma::mat b){
//   arma::mat dm = est_B_m - oldest_B_m;
//   int nrb = B.n_rows;
//   arma::mat invW;
//   arma::mat Z;
//   while((as_scalar(trans(dm) * dm) / nrb) > pow(tol_pql, 2) &&
//   (iteration_m <= maxit_pql)){
//     iteration_m = iteration_m + 1;
//     oldest_B_m = est_B_m;
//     arma::vec pq = 1 / (mu % (1 - mu));
//     invW = mat(diagmat(pq));
//     arma::mat V = invW + C;
//     if(V.has_inf() || (rcond(V) < pow(static_cast<double>(10), static_cast<double>(-10)))){
//       // Rcout << rcondflag << " " ;
//       rcondflag = rcondflag + 1;
//       B = mat(size(B_init));
//       B.fill(0.001);
//       b = mat(n, 1);
//       b.fill(0.0);
//       arma::mat beta = join_cols(B, b);
//       mu = exp(X * B) / (1 + exp(X * B));
//       arma::mat oldest_B_m = mat(est_B.n_rows, 1);
//       oldest_B_m.fill(1000000);
//       arma::vec pq2 = mu % (1 - mu);
//       invW.diag() = pq2;
//       V = invW + C;
//     }
//     arma::mat invV = inv(V);
//     Z = X * B + b + (y - mu)/(mu % (1 - mu));
//     arma::mat denom = trans(X) * invV * X;
//     arma::mat num = trans(X) * invV * Z;
//     B = solve(denom, num);
//     b = C * invV * (Z - X * B);
//     arma::mat beta2 = join_cols(B, b);
//     mu = exp(XX * beta2) / (1 + exp(XX * beta2));
//     est_B_m = mat(B);
//   }
//   
//   List out = List::create(
//     _["Z"] = Z, _["B"] = B, _["b"] = b, _["mu"] = mu, _["invW"] = invW,
//     _["est.B.m"] = est_B_m, _["rcondflag"] = rcondflag
//   );
//   
//   return out;
// }
// 
// // [[Rcpp::export]]
// List binpglmm_inter_while_cpp2(
//     arma::mat est_B_m, 
//     arma::mat B, 
//     arma::vec mu, 
//     const arma::mat& C, 
//     int rcondflag, 
//     const arma::mat& B_init, 
//     const arma::mat& X, 
//     const arma::mat& XX,
//     const arma::mat& est_B, 
//     const arma::vec& y, 
//     int n, arma::mat b){
//     arma::vec pq = 1 / (mu % (1 - mu));
//     arma::mat invW = mat(diagmat(pq));
//     arma::mat V = invW + C;
//     if(V.has_inf() || (rcond(V) < pow(static_cast<double>(10), static_cast<double>(-10)))){
//       // Rcout << rcondflag << " " ;
//       rcondflag = rcondflag + 1;
//       B = mat(size(B_init));
//       B.fill(0.001);
//       b = mat(n, 1);
//       b.fill(0.0);
//       arma::mat beta = join_cols(B, b);
//       mu = exp(X * B) / (1 + exp(X * B));
//       arma::mat oldest_B_m = mat(est_B.n_rows, 1);
//       oldest_B_m.fill(1000000);
//       arma::vec pq2 = mu % (1 - mu);
//       invW.diag() = pq2;
//       V = invW + C;
//     }
//     arma::mat invV = inv(V);
//     arma::mat Z = X * B + b + (y - mu)/(mu % (1 - mu));
//     arma::mat denom = trans(X) * invV * X;
//     arma::mat num = trans(X) * invV * Z;
//     B = solve(denom, num);
//     b = C * invV * (Z - X * B);
//     arma::mat beta2 = join_cols(B, b);
//     mu = exp(XX * beta2) / (1 + exp(XX * beta2));
//     est_B_m = mat(B);
//   // }
//   
//   List out = List::create(
//     _["Z"] = Z, _["B"] = B, _["b"] = b, _["mu"] = mu, _["invW"] = invW,
//       _["est.B.m"] = est_B_m, _["rcondflag"] = rcondflag
//   );
//   
//   return out;
// }
// 
// 
// /*** R
// # pglmm.reml(par = s2, tinvW = invW, tH = H, tVphy = as.matrix(Vphy), tX = X)
// # pglmm_reml_cpp(par = s2, tinvW = as.matrix(invW), tH = H, tVphy = as.matrix(Vphy), tX = X)
// # xx1 = binpglmm_inter_while_cpp(est.B.m, oldest.B.m, B, tol.pql, iteration.m, maxit.pql,
// #                          mu, C, rcondflag, B.init, X, XX, est.B, y, n, b)
// # xx2 = binpglmm_inter_while(est.B.m, oldest.B.m, B, tol.pql, iteration.m, maxit.pql, mu,
// #                      C, rcondflag, B.init, X, XX, est.B, y, n, b)
// # # testthat::expect_equivalent(xx1, xx2)
// # round(xx1$Z[,1], 5) == round(unname(xx2$Z[,1]), 5)
// # xx1$B
// # xx2$B
// */
