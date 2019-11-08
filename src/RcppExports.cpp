// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// pglmm_reml_cpp
double pglmm_reml_cpp(arma::vec par, const arma::mat& tinvW, const arma::mat& tH, const arma::mat& tVphy, const arma::mat& tX);
RcppExport SEXP _phyr_pglmm_reml_cpp(SEXP parSEXP, SEXP tinvWSEXP, SEXP tHSEXP, SEXP tVphySEXP, SEXP tXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tinvW(tinvWSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tH(tHSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tVphy(tVphySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tX(tXSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_reml_cpp(par, tinvW, tH, tVphy, tX));
    return rcpp_result_gen;
END_RCPP
}
// binpglmm_inter_while_cpp
List binpglmm_inter_while_cpp(arma::mat est_B_m, arma::mat oldest_B_m, arma::mat B, const double& tol_pql, int iteration_m, const int& maxit_pql, arma::vec mu, const arma::mat& C, int rcondflag, const arma::mat& B_init, const arma::mat& X, const arma::mat& XX, const arma::mat& est_B, const arma::vec& y, int n, arma::mat b);
RcppExport SEXP _phyr_binpglmm_inter_while_cpp(SEXP est_B_mSEXP, SEXP oldest_B_mSEXP, SEXP BSEXP, SEXP tol_pqlSEXP, SEXP iteration_mSEXP, SEXP maxit_pqlSEXP, SEXP muSEXP, SEXP CSEXP, SEXP rcondflagSEXP, SEXP B_initSEXP, SEXP XSEXP, SEXP XXSEXP, SEXP est_BSEXP, SEXP ySEXP, SEXP nSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type est_B_m(est_B_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type oldest_B_m(oldest_B_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol_pql(tol_pqlSEXP);
    Rcpp::traits::input_parameter< int >::type iteration_m(iteration_mSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxit_pql(maxit_pqlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type rcondflag(rcondflagSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_init(B_initSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type est_B(est_BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(binpglmm_inter_while_cpp(est_B_m, oldest_B_m, B, tol_pql, iteration_m, maxit_pql, mu, C, rcondflag, B_init, X, XX, est_B, y, n, b));
    return rcpp_result_gen;
END_RCPP
}
// binpglmm_inter_while_cpp2
List binpglmm_inter_while_cpp2(arma::mat est_B_m, arma::mat B, arma::vec mu, const arma::mat& C, int rcondflag, const arma::mat& B_init, const arma::mat& X, const arma::mat& XX, const arma::mat& est_B, const arma::vec& y, int n, arma::mat b);
RcppExport SEXP _phyr_binpglmm_inter_while_cpp2(SEXP est_B_mSEXP, SEXP BSEXP, SEXP muSEXP, SEXP CSEXP, SEXP rcondflagSEXP, SEXP B_initSEXP, SEXP XSEXP, SEXP XXSEXP, SEXP est_BSEXP, SEXP ySEXP, SEXP nSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type est_B_m(est_B_mSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type rcondflag(rcondflagSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B_init(B_initSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type est_B(est_BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(binpglmm_inter_while_cpp2(est_B_m, B, mu, C, rcondflag, B_init, X, XX, est_B, y, n, b));
    return rcpp_result_gen;
END_RCPP
}
// cor_phylo_LL
double cor_phylo_LL(NumericVector par, SEXP xptr);
RcppExport SEXP _phyr_cor_phylo_LL(SEXP parSEXP, SEXP xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_phylo_LL(par, xptr));
    return rcpp_result_gen;
END_RCPP
}
// cor_phylo_cpp
List cor_phylo_cpp(const arma::mat& X, const std::vector<arma::mat>& U, const arma::mat& M, const arma::mat& Vphy_, const bool& REML, const bool& constrain_d, const double& lower_d, const bool& verbose, const double& rcond_threshold, const double& rel_tol, const int& max_iter, const std::string& method, const bool& no_corr, const uint_fast32_t& boot, const std::string& keep_boots, const std::vector<double>& sann);
RcppExport SEXP _phyr_cor_phylo_cpp(SEXP XSEXP, SEXP USEXP, SEXP MSEXP, SEXP Vphy_SEXP, SEXP REMLSEXP, SEXP constrain_dSEXP, SEXP lower_dSEXP, SEXP verboseSEXP, SEXP rcond_thresholdSEXP, SEXP rel_tolSEXP, SEXP max_iterSEXP, SEXP methodSEXP, SEXP no_corrSEXP, SEXP bootSEXP, SEXP keep_bootsSEXP, SEXP sannSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Vphy_(Vphy_SEXP);
    Rcpp::traits::input_parameter< const bool& >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< const bool& >::type constrain_d(constrain_dSEXP);
    Rcpp::traits::input_parameter< const double& >::type lower_d(lower_dSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const double& >::type rcond_threshold(rcond_thresholdSEXP);
    Rcpp::traits::input_parameter< const double& >::type rel_tol(rel_tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const bool& >::type no_corr(no_corrSEXP);
    Rcpp::traits::input_parameter< const uint_fast32_t& >::type boot(bootSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type keep_boots(keep_bootsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type sann(sannSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_phylo_cpp(X, U, M, Vphy_, REML, constrain_d, lower_d, verbose, rcond_threshold, rel_tol, max_iter, method, no_corr, boot, keep_boots, sann));
    return rcpp_result_gen;
END_RCPP
}
// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP _phyr_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}
// predict_cpp
NumericVector predict_cpp(int n, const arma::vec& nsr, int reps, const arma::mat& V);
RcppExport SEXP _phyr_predict_cpp(SEXP nSEXP, SEXP nsrSEXP, SEXP repsSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nsr(nsrSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_cpp(n, nsr, reps, V));
    return rcpp_result_gen;
END_RCPP
}
// pcd2_loop
List pcd2_loop(arma::vec SSii, arma::vec nsr, double SCii, const arma::mat& comm, const arma::mat& V, int nsp_pool, bool verbose);
RcppExport SEXP _phyr_pcd2_loop(SEXP SSiiSEXP, SEXP nsrSEXP, SEXP SCiiSEXP, SEXP commSEXP, SEXP VSEXP, SEXP nsp_poolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type SSii(SSiiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nsr(nsrSEXP);
    Rcpp::traits::input_parameter< double >::type SCii(SCiiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type comm(commSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type nsp_pool(nsp_poolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pcd2_loop(SSii, nsr, SCii, comm, V, nsp_pool, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_iV_logdetV_cpp
List pglmm_iV_logdetV_cpp(NumericVector par, arma::vec mu, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, bool logdet, const std::string family, arma::vec totalSize);
RcppExport SEXP _phyr_pglmm_iV_logdetV_cpp(SEXP parSEXP, SEXP muSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP nestedSEXP, SEXP logdetSEXP, SEXP familySEXP, SEXP totalSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type logdet(logdetSEXP);
    Rcpp::traits::input_parameter< const std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type totalSize(totalSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_iV_logdetV_cpp(par, mu, Zt, St, nested, logdet, family, totalSize));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_V
arma::sp_mat pglmm_V(NumericVector par, const arma::sp_mat& Zt, const arma::sp_mat& St, arma::vec mu, const List& nested, bool missing_mu, const std::string family, arma::vec totalSize);
RcppExport SEXP _phyr_pglmm_V(SEXP parSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP muSEXP, SEXP nestedSEXP, SEXP missing_muSEXP, SEXP familySEXP, SEXP totalSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type missing_mu(missing_muSEXP);
    Rcpp::traits::input_parameter< const std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type totalSize(totalSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_V(par, Zt, St, mu, nested, missing_mu, family, totalSize));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_LL_cpp
double pglmm_LL_cpp(NumericVector par, const arma::vec& H, const arma::mat& X, const arma::sp_mat& Zt, const arma::sp_mat& St, const arma::vec& mu, const List& nested, bool REML, bool verbose, const std::string family, arma::vec totalSize);
RcppExport SEXP _phyr_pglmm_LL_cpp(SEXP parSEXP, SEXP HSEXP, SEXP XSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP muSEXP, SEXP nestedSEXP, SEXP REMLSEXP, SEXP verboseSEXP, SEXP familySEXP, SEXP totalSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type totalSize(totalSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_LL_cpp(par, H, X, Zt, St, mu, nested, REML, verbose, family, totalSize));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_internal_cpp
List pglmm_internal_cpp(const arma::mat& X, const arma::vec& Y, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, const bool REML, const bool verbose, const int n, const int p, const int q, const int maxit, const double reltol, const double tol_pql, const double maxit_pql, const std::string optimizer, arma::mat B_init, arma::vec ss, const std::string family, arma::vec totalSize);
RcppExport SEXP _phyr_pglmm_internal_cpp(SEXP XSEXP, SEXP YSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP nestedSEXP, SEXP REMLSEXP, SEXP verboseSEXP, SEXP nSEXP, SEXP pSEXP, SEXP qSEXP, SEXP maxitSEXP, SEXP reltolSEXP, SEXP tol_pqlSEXP, SEXP maxit_pqlSEXP, SEXP optimizerSEXP, SEXP B_initSEXP, SEXP ssSEXP, SEXP familySEXP, SEXP totalSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< const bool >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_pql(tol_pqlSEXP);
    Rcpp::traits::input_parameter< const double >::type maxit_pql(maxit_pqlSEXP);
    Rcpp::traits::input_parameter< const std::string >::type optimizer(optimizerSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B_init(B_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< const std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type totalSize(totalSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_internal_cpp(X, Y, Zt, St, nested, REML, verbose, n, p, q, maxit, reltol, tol_pql, maxit_pql, optimizer, B_init, ss, family, totalSize));
    return rcpp_result_gen;
END_RCPP
}
// sexp_type
int sexp_type(SEXP x);
RcppExport SEXP _phyr_sexp_type(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sexp_type(x));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_gaussian_predict
arma::vec pglmm_gaussian_predict(const arma::mat& iV, const arma::mat& H);
RcppExport SEXP _phyr_pglmm_gaussian_predict(SEXP iVSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type iV(iVSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_gaussian_predict(iV, H));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_gaussian_LL_cpp
double pglmm_gaussian_LL_cpp(NumericVector par, const arma::mat& X, const arma::vec& Y, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, bool REML, bool verbose);
RcppExport SEXP _phyr_pglmm_gaussian_LL_cpp(SEXP parSEXP, SEXP XSEXP, SEXP YSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP nestedSEXP, SEXP REMLSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_gaussian_LL_cpp(par, X, Y, Zt, St, nested, REML, verbose));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_gaussian_LL_calc_cpp
List pglmm_gaussian_LL_calc_cpp(NumericVector par, const arma::mat& X, const arma::vec& Y, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, bool REML);
RcppExport SEXP _phyr_pglmm_gaussian_LL_calc_cpp(SEXP parSEXP, SEXP XSEXP, SEXP YSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP nestedSEXP, SEXP REMLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_gaussian_LL_calc_cpp(par, X, Y, Zt, St, nested, REML));
    return rcpp_result_gen;
END_RCPP
}
// pglmm_gaussian_internal_cpp
Rcpp::List pglmm_gaussian_internal_cpp(NumericVector par, const arma::mat& X, const arma::vec& Y, const arma::sp_mat& Zt, const arma::sp_mat& St, const List& nested, bool REML, bool verbose, std::string optimizer, int maxit, double reltol, int q, int n, int p, const double Pi);
RcppExport SEXP _phyr_pglmm_gaussian_internal_cpp(SEXP parSEXP, SEXP XSEXP, SEXP YSEXP, SEXP ZtSEXP, SEXP StSEXP, SEXP nestedSEXP, SEXP REMLSEXP, SEXP verboseSEXP, SEXP optimizerSEXP, SEXP maxitSEXP, SEXP reltolSEXP, SEXP qSEXP, SEXP nSEXP, SEXP pSEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Zt(ZtSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type St(StSEXP);
    Rcpp::traits::input_parameter< const List& >::type nested(nestedSEXP);
    Rcpp::traits::input_parameter< bool >::type REML(REMLSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< std::string >::type optimizer(optimizerSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(pglmm_gaussian_internal_cpp(par, X, Y, Zt, St, nested, REML, verbose, optimizer, maxit, reltol, q, n, p, Pi));
    return rcpp_result_gen;
END_RCPP
}
// which2
IntegerVector which2(const LogicalVector x);
RcppExport SEXP _phyr_which2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which2(x));
    return rcpp_result_gen;
END_RCPP
}
// vcv_loop
arma::mat vcv_loop(NumericVector& xx, const int& n, const IntegerVector& e1, const IntegerVector& e2, const NumericVector& EL, const List& pp, bool corr);
RcppExport SEXP _phyr_vcv_loop(SEXP xxSEXP, SEXP nSEXP, SEXP e1SEXP, SEXP e2SEXP, SEXP ELSEXP, SEXP ppSEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type e1(e1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type e2(e2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type EL(ELSEXP);
    Rcpp::traits::input_parameter< const List& >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< bool >::type corr(corrSEXP);
    rcpp_result_gen = Rcpp::wrap(vcv_loop(xx, n, e1, e2, EL, pp, corr));
    return rcpp_result_gen;
END_RCPP
}
// cov2cor_cpp
void cov2cor_cpp(arma::mat& vcv);
RcppExport SEXP _phyr_cov2cor_cpp(SEXP vcvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type vcv(vcvSEXP);
    cov2cor_cpp(vcv);
    return R_NilValue;
END_RCPP
}
// pse_cpp
NumericVector pse_cpp(const NumericMatrix& comm, const arma::mat& Cmatrix);
RcppExport SEXP _phyr_pse_cpp(SEXP commSEXP, SEXP CmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type comm(commSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Cmatrix(CmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(pse_cpp(comm, Cmatrix));
    return rcpp_result_gen;
END_RCPP
}
// psv_cpp
DataFrame psv_cpp(const NumericMatrix& comm, const arma::mat& Cmatrix, const bool compute_var);
RcppExport SEXP _phyr_psv_cpp(SEXP commSEXP, SEXP CmatrixSEXP, SEXP compute_varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type comm(commSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Cmatrix(CmatrixSEXP);
    Rcpp::traits::input_parameter< const bool >::type compute_var(compute_varSEXP);
    rcpp_result_gen = Rcpp::wrap(psv_cpp(comm, Cmatrix, compute_var));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_phyr_pglmm_reml_cpp", (DL_FUNC) &_phyr_pglmm_reml_cpp, 5},
    {"_phyr_binpglmm_inter_while_cpp", (DL_FUNC) &_phyr_binpglmm_inter_while_cpp, 16},
    {"_phyr_binpglmm_inter_while_cpp2", (DL_FUNC) &_phyr_binpglmm_inter_while_cpp2, 12},
    {"_phyr_cor_phylo_LL", (DL_FUNC) &_phyr_cor_phylo_LL, 2},
    {"_phyr_cor_phylo_cpp", (DL_FUNC) &_phyr_cor_phylo_cpp, 16},
    {"_phyr_set_seed", (DL_FUNC) &_phyr_set_seed, 1},
    {"_phyr_predict_cpp", (DL_FUNC) &_phyr_predict_cpp, 4},
    {"_phyr_pcd2_loop", (DL_FUNC) &_phyr_pcd2_loop, 7},
    {"_phyr_pglmm_iV_logdetV_cpp", (DL_FUNC) &_phyr_pglmm_iV_logdetV_cpp, 8},
    {"_phyr_pglmm_V", (DL_FUNC) &_phyr_pglmm_V, 8},
    {"_phyr_pglmm_LL_cpp", (DL_FUNC) &_phyr_pglmm_LL_cpp, 11},
    {"_phyr_pglmm_internal_cpp", (DL_FUNC) &_phyr_pglmm_internal_cpp, 19},
    {"_phyr_sexp_type", (DL_FUNC) &_phyr_sexp_type, 1},
    {"_phyr_pglmm_gaussian_predict", (DL_FUNC) &_phyr_pglmm_gaussian_predict, 2},
    {"_phyr_pglmm_gaussian_LL_cpp", (DL_FUNC) &_phyr_pglmm_gaussian_LL_cpp, 8},
    {"_phyr_pglmm_gaussian_LL_calc_cpp", (DL_FUNC) &_phyr_pglmm_gaussian_LL_calc_cpp, 7},
    {"_phyr_pglmm_gaussian_internal_cpp", (DL_FUNC) &_phyr_pglmm_gaussian_internal_cpp, 15},
    {"_phyr_which2", (DL_FUNC) &_phyr_which2, 1},
    {"_phyr_vcv_loop", (DL_FUNC) &_phyr_vcv_loop, 7},
    {"_phyr_cov2cor_cpp", (DL_FUNC) &_phyr_cov2cor_cpp, 1},
    {"_phyr_pse_cpp", (DL_FUNC) &_phyr_pse_cpp, 2},
    {"_phyr_psv_cpp", (DL_FUNC) &_phyr_psv_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_phyr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
