// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// h_tilde
NumericVector h_tilde(NumericVector X1, NumericVector X2, NumericVector X3, NumericVector X4, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGGM_h_tilde(SEXP X1SEXP, SEXP X2SEXP, SEXP X3SEXP, SEXP X4SEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X3(X3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X4(X4SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(h_tilde(X1, X2, X3, X4, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// h
NumericVector h(List L, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2, IntegerMatrix perm);
RcppExport SEXP _TestGGM_h(SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(h(L, ind_eq, ind_ineq1, ind_ineq2, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_H
NumericMatrix calculate_H(NumericMatrix X, IntegerMatrix indices_U, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGGM_calculate_H(SEXP XSEXP, SEXP indices_USEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indices_U(indices_USEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// calculate_H_not_symmetric
NumericMatrix calculate_H_not_symmetric(NumericMatrix X, IntegerMatrix indices_U, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGGM_calculate_H_not_symmetric(SEXP XSEXP, SEXP indices_USEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indices_U(indices_USEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_H_not_symmetric(X, indices_U, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// calculate_g
NumericVector calculate_g(NumericMatrix X, int i, int L, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2, IntegerMatrix perm);
RcppExport SEXP _TestGGM_calculate_g(SEXP XSEXP, SEXP iSEXP, SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_g(X, i, L, ind_eq, ind_ineq1, ind_ineq2, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G
NumericMatrix calculate_G(NumericMatrix X, int L, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGGM_calculate_G(SEXP XSEXP, SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_G(X, L, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// compute_S_small
IntegerMatrix compute_S_small(int n, int i, int L, int K);
RcppExport SEXP _TestGGM_compute_S_small(SEXP nSEXP, SEXP iSEXP, SEXP LSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_S_small(n, i, L, K));
    return rcpp_result_gen;
END_RCPP
}
// g_small
NumericVector g_small(NumericMatrix X, int i, int L, int n1, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2, IntegerMatrix perm);
RcppExport SEXP _TestGGM_g_small(SEXP XSEXP, SEXP iSEXP, SEXP LSEXP, SEXP n1SEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(g_small(X, i, L, n1, ind_eq, ind_ineq1, ind_ineq2, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G_small
NumericMatrix calculate_G_small(NumericMatrix X, int L, int n1, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGGM_calculate_G_small(SEXP XSEXP, SEXP LSEXP, SEXP n1SEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_G_small(X, L, n1, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// h_tilde_eq
NumericVector h_tilde_eq(NumericVector X1, NumericVector X2, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_h_tilde_eq(SEXP X1SEXP, SEXP X2SEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(h_tilde_eq(X1, X2, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// h_eq
NumericVector h_eq(List L, IntegerMatrix ind_eq, IntegerMatrix perm);
RcppExport SEXP _TestGGM_h_eq(SEXP LSEXP, SEXP ind_eqSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(h_eq(L, ind_eq, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_H_eq
NumericMatrix calculate_H_eq(NumericMatrix X, IntegerMatrix indices_U, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_calculate_H_eq(SEXP XSEXP, SEXP indices_USEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indices_U(indices_USEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_H_eq(X, indices_U, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// calculate_H_not_symmetric_eq
NumericMatrix calculate_H_not_symmetric_eq(NumericMatrix X, IntegerMatrix indices_U, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_calculate_H_not_symmetric_eq(SEXP XSEXP, SEXP indices_USEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indices_U(indices_USEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_H_not_symmetric_eq(X, indices_U, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// g_eq
NumericVector g_eq(NumericMatrix X, int i, int L, IntegerMatrix ind_eq, IntegerMatrix perm);
RcppExport SEXP _TestGGM_g_eq(SEXP XSEXP, SEXP iSEXP, SEXP LSEXP, SEXP ind_eqSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(g_eq(X, i, L, ind_eq, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G_eq
NumericMatrix calculate_G_eq(NumericMatrix X, int L, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_calculate_G_eq(SEXP XSEXP, SEXP LSEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_G_eq(X, L, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap
NumericMatrix bootstrap(int E, NumericMatrix H_centered);
RcppExport SEXP _TestGGM_bootstrap(SEXP ESEXP, SEXP H_centeredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H_centered(H_centeredSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap(E, H_centered));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_m_dep
NumericMatrix bootstrap_m_dep(int E, int B, int omega, NumericMatrix Y_centered);
RcppExport SEXP _TestGGM_bootstrap_m_dep(SEXP ESEXP, SEXP BSEXP, SEXP omegaSEXP, SEXP Y_centeredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_centered(Y_centeredSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_m_dep(E, B, omega, Y_centered));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_U
List bootstrap_U(int E, int r, NumericMatrix H_centered, NumericMatrix G_centered);
RcppExport SEXP _TestGGM_bootstrap_U(SEXP ESEXP, SEXP rSEXP, SEXP H_centeredSEXP, SEXP G_centeredSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H_centered(H_centeredSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G_centered(G_centeredSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_U(E, r, H_centered, G_centered));
    return rcpp_result_gen;
END_RCPP
}
// expectation
double expectation(NumericMatrix S, IntegerVector j);
RcppExport SEXP _TestGGM_expectation(SEXP SSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(expectation(S, j));
    return rcpp_result_gen;
END_RCPP
}
// fourth_mom
double fourth_mom(NumericMatrix S, int a, int b, int c, int d);
RcppExport SEXP _TestGGM_fourth_mom(SEXP SSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(fourth_mom(S, a, b, c, d));
    return rcpp_result_gen;
END_RCPP
}
// cov_grouping
NumericMatrix cov_grouping(NumericMatrix S, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_cov_grouping(SEXP SSEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_grouping(S, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// cov_run_over
NumericMatrix cov_run_over(NumericMatrix S, IntegerMatrix ind_eq);
RcppExport SEXP _TestGGM_cov_run_over(SEXP SSEXP, SEXP ind_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_run_over(S, ind_eq));
    return rcpp_result_gen;
END_RCPP
}
// h_tilde_fac
NumericVector h_tilde_fac(NumericVector X1, NumericVector X2, NumericVector X3, IntegerMatrix ind);
RcppExport SEXP _TestGGM_h_tilde_fac(SEXP X1SEXP, SEXP X2SEXP, SEXP X3SEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type X3(X3SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(h_tilde_fac(X1, X2, X3, ind));
    return rcpp_result_gen;
END_RCPP
}
// h_fac
NumericVector h_fac(List L, IntegerMatrix ind_minors, IntegerMatrix perm);
RcppExport SEXP _TestGGM_h_fac(SEXP LSEXP, SEXP ind_minorsSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_minors(ind_minorsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(h_fac(L, ind_minors, perm));
    return rcpp_result_gen;
END_RCPP
}
// H_factors
NumericMatrix H_factors(NumericMatrix X, IntegerMatrix indices, IntegerMatrix ind_minors);
RcppExport SEXP _TestGGM_H_factors(SEXP XSEXP, SEXP indicesSEXP, SEXP ind_minorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_minors(ind_minorsSEXP);
    rcpp_result_gen = Rcpp::wrap(H_factors(X, indices, ind_minors));
    return rcpp_result_gen;
END_RCPP
}
// g_fac
NumericVector g_fac(NumericMatrix X, int i, int L, IntegerMatrix ind, IntegerMatrix perm);
RcppExport SEXP _TestGGM_g_fac(SEXP XSEXP, SEXP iSEXP, SEXP LSEXP, SEXP indSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind(indSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(g_fac(X, i, L, ind, perm));
    return rcpp_result_gen;
END_RCPP
}
// G_factors
NumericMatrix G_factors(NumericMatrix X, IntegerVector S1, int L, IntegerMatrix ind);
RcppExport SEXP _TestGGM_G_factors(SEXP XSEXP, SEXP S1SEXP, SEXP LSEXP, SEXP indSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind(indSEXP);
    rcpp_result_gen = Rcpp::wrap(G_factors(X, S1, L, ind));
    return rcpp_result_gen;
END_RCPP
}
// random_combs
List random_combs(int n, int k, int nr);
RcppExport SEXP _TestGGM_random_combs(SEXP nSEXP, SEXP kSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(random_combs(n, k, nr));
    return rcpp_result_gen;
END_RCPP
}
// findn
int findn(int N, int D);
RcppExport SEXP _TestGGM_findn(SEXP NSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(findn(N, D));
    return rcpp_result_gen;
END_RCPP
}
// compute_S
IntegerMatrix compute_S(int n, int i, int L);
RcppExport SEXP _TestGGM_compute_S(SEXP nSEXP, SEXP iSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_S(n, i, L));
    return rcpp_result_gen;
END_RCPP
}
// permutations
IntegerMatrix permutations(int n);
RcppExport SEXP _TestGGM_permutations(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(permutations(n));
    return rcpp_result_gen;
END_RCPP
}
// update_param
List update_param(NumericMatrix S, NumericMatrix edges, int nr_obs);
RcppExport SEXP _TestGGM_update_param(SEXP SSEXP, SEXP edgesSEXP, SEXP nr_obsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< int >::type nr_obs(nr_obsSEXP);
    rcpp_result_gen = Rcpp::wrap(update_param(S, edges, nr_obs));
    return rcpp_result_gen;
END_RCPP
}
// cov_from_graph_large
NumericMatrix cov_from_graph_large(NumericVector Omega, NumericVector Rho, List paths);
RcppExport SEXP _TestGGM_cov_from_graph_large(SEXP OmegaSEXP, SEXP RhoSEXP, SEXP pathsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< List >::type paths(pathsSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_from_graph_large(Omega, Rho, paths));
    return rcpp_result_gen;
END_RCPP
}
// binomialCoefficient
int binomialCoefficient(int n, int k);
RcppExport SEXP _TestGGM_binomialCoefficient(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(binomialCoefficient(n, k));
    return rcpp_result_gen;
END_RCPP
}
// two_subsets
IntegerMatrix two_subsets(int n);
RcppExport SEXP _TestGGM_two_subsets(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(two_subsets(n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TestGGM_h_tilde", (DL_FUNC) &_TestGGM_h_tilde, 7},
    {"_TestGGM_h", (DL_FUNC) &_TestGGM_h, 5},
    {"_TestGGM_calculate_H", (DL_FUNC) &_TestGGM_calculate_H, 5},
    {"_TestGGM_calculate_H_not_symmetric", (DL_FUNC) &_TestGGM_calculate_H_not_symmetric, 5},
    {"_TestGGM_calculate_g", (DL_FUNC) &_TestGGM_calculate_g, 7},
    {"_TestGGM_calculate_G", (DL_FUNC) &_TestGGM_calculate_G, 5},
    {"_TestGGM_compute_S_small", (DL_FUNC) &_TestGGM_compute_S_small, 4},
    {"_TestGGM_g_small", (DL_FUNC) &_TestGGM_g_small, 8},
    {"_TestGGM_calculate_G_small", (DL_FUNC) &_TestGGM_calculate_G_small, 6},
    {"_TestGGM_h_tilde_eq", (DL_FUNC) &_TestGGM_h_tilde_eq, 3},
    {"_TestGGM_h_eq", (DL_FUNC) &_TestGGM_h_eq, 3},
    {"_TestGGM_calculate_H_eq", (DL_FUNC) &_TestGGM_calculate_H_eq, 3},
    {"_TestGGM_calculate_H_not_symmetric_eq", (DL_FUNC) &_TestGGM_calculate_H_not_symmetric_eq, 3},
    {"_TestGGM_g_eq", (DL_FUNC) &_TestGGM_g_eq, 5},
    {"_TestGGM_calculate_G_eq", (DL_FUNC) &_TestGGM_calculate_G_eq, 3},
    {"_TestGGM_bootstrap", (DL_FUNC) &_TestGGM_bootstrap, 2},
    {"_TestGGM_bootstrap_m_dep", (DL_FUNC) &_TestGGM_bootstrap_m_dep, 4},
    {"_TestGGM_bootstrap_U", (DL_FUNC) &_TestGGM_bootstrap_U, 4},
    {"_TestGGM_expectation", (DL_FUNC) &_TestGGM_expectation, 2},
    {"_TestGGM_fourth_mom", (DL_FUNC) &_TestGGM_fourth_mom, 5},
    {"_TestGGM_cov_grouping", (DL_FUNC) &_TestGGM_cov_grouping, 2},
    {"_TestGGM_cov_run_over", (DL_FUNC) &_TestGGM_cov_run_over, 2},
    {"_TestGGM_h_tilde_fac", (DL_FUNC) &_TestGGM_h_tilde_fac, 4},
    {"_TestGGM_h_fac", (DL_FUNC) &_TestGGM_h_fac, 3},
    {"_TestGGM_H_factors", (DL_FUNC) &_TestGGM_H_factors, 3},
    {"_TestGGM_g_fac", (DL_FUNC) &_TestGGM_g_fac, 5},
    {"_TestGGM_G_factors", (DL_FUNC) &_TestGGM_G_factors, 4},
    {"_TestGGM_random_combs", (DL_FUNC) &_TestGGM_random_combs, 3},
    {"_TestGGM_findn", (DL_FUNC) &_TestGGM_findn, 2},
    {"_TestGGM_compute_S", (DL_FUNC) &_TestGGM_compute_S, 3},
    {"_TestGGM_permutations", (DL_FUNC) &_TestGGM_permutations, 1},
    {"_TestGGM_update_param", (DL_FUNC) &_TestGGM_update_param, 3},
    {"_TestGGM_cov_from_graph_large", (DL_FUNC) &_TestGGM_cov_from_graph_large, 3},
    {"_TestGGM_binomialCoefficient", (DL_FUNC) &_TestGGM_binomialCoefficient, 2},
    {"_TestGGM_two_subsets", (DL_FUNC) &_TestGGM_two_subsets, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_TestGGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
