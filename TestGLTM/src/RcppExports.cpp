// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_S
IntegerMatrix compute_S(int n, int i, int L);
RcppExport SEXP _TestGLTM_compute_S(SEXP nSEXP, SEXP iSEXP, SEXP LSEXP) {
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
RcppExport SEXP _TestGLTM_permutations(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(permutations(n));
    return rcpp_result_gen;
END_RCPP
}
// h_tilde
NumericVector h_tilde(NumericVector X1, NumericVector X2, NumericVector X3, NumericVector X4, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGLTM_h_tilde(SEXP X1SEXP, SEXP X2SEXP, SEXP X3SEXP, SEXP X4SEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
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
RcppExport SEXP _TestGLTM_h(SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
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
RcppExport SEXP _TestGLTM_calculate_H(SEXP XSEXP, SEXP indices_USEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
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
// g
NumericVector g(NumericMatrix X, int i, int L, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2, IntegerMatrix perm);
RcppExport SEXP _TestGLTM_g(SEXP XSEXP, SEXP iSEXP, SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(g(X, i, L, ind_eq, ind_ineq1, ind_ineq2, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_G
NumericMatrix calculate_G(NumericMatrix X, int L, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGLTM_calculate_G(SEXP XSEXP, SEXP LSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
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
// findn
int findn(int N, int D);
RcppExport SEXP _TestGLTM_findn(SEXP NSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(findn(N, D));
    return rcpp_result_gen;
END_RCPP
}
// calculate_Y_independent
NumericMatrix calculate_Y_independent(NumericMatrix X, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGLTM_calculate_Y_independent(SEXP XSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_Y_independent(X, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// calculate_Y_symmetric
NumericMatrix calculate_Y_symmetric(NumericMatrix X, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2, IntegerMatrix perm);
RcppExport SEXP _TestGLTM_calculate_Y_symmetric(SEXP XSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_Y_symmetric(X, ind_eq, ind_ineq1, ind_ineq2, perm));
    return rcpp_result_gen;
END_RCPP
}
// calculate_Y_m_dep
NumericMatrix calculate_Y_m_dep(NumericMatrix X, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2);
RcppExport SEXP _TestGLTM_calculate_Y_m_dep(SEXP XSEXP, SEXP ind_eqSEXP, SEXP ind_ineq1SEXP, SEXP ind_ineq2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_eq(ind_eqSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq1(ind_ineq1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type ind_ineq2(ind_ineq2SEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_Y_m_dep(X, ind_eq, ind_ineq1, ind_ineq2));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_independent
NumericVector bootstrap_independent(int E, NumericVector standardizer, NumericMatrix Y_centered, int p_eq);
RcppExport SEXP _TestGLTM_bootstrap_independent(SEXP ESEXP, SEXP standardizerSEXP, SEXP Y_centeredSEXP, SEXP p_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type standardizer(standardizerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_centered(Y_centeredSEXP);
    Rcpp::traits::input_parameter< int >::type p_eq(p_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_independent(E, standardizer, Y_centered, p_eq));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_m_dep
NumericVector bootstrap_m_dep(int E, int B, int omega, NumericVector standardizer, NumericMatrix Y_centered, int p_eq);
RcppExport SEXP _TestGLTM_bootstrap_m_dep(SEXP ESEXP, SEXP BSEXP, SEXP omegaSEXP, SEXP standardizerSEXP, SEXP Y_centeredSEXP, SEXP p_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type standardizer(standardizerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y_centered(Y_centeredSEXP);
    Rcpp::traits::input_parameter< int >::type p_eq(p_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_m_dep(E, B, omega, standardizer, Y_centered, p_eq));
    return rcpp_result_gen;
END_RCPP
}
// bootstrap_U
List bootstrap_U(int E, NumericMatrix H_centered, NumericMatrix G_centered, int N);
RcppExport SEXP _TestGLTM_bootstrap_U(SEXP ESEXP, SEXP H_centeredSEXP, SEXP G_centeredSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type H_centered(H_centeredSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G_centered(G_centeredSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrap_U(E, H_centered, G_centered, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TestGLTM_compute_S", (DL_FUNC) &_TestGLTM_compute_S, 3},
    {"_TestGLTM_permutations", (DL_FUNC) &_TestGLTM_permutations, 1},
    {"_TestGLTM_h_tilde", (DL_FUNC) &_TestGLTM_h_tilde, 7},
    {"_TestGLTM_h", (DL_FUNC) &_TestGLTM_h, 5},
    {"_TestGLTM_calculate_H", (DL_FUNC) &_TestGLTM_calculate_H, 5},
    {"_TestGLTM_g", (DL_FUNC) &_TestGLTM_g, 7},
    {"_TestGLTM_calculate_G", (DL_FUNC) &_TestGLTM_calculate_G, 5},
    {"_TestGLTM_findn", (DL_FUNC) &_TestGLTM_findn, 2},
    {"_TestGLTM_calculate_Y_independent", (DL_FUNC) &_TestGLTM_calculate_Y_independent, 4},
    {"_TestGLTM_calculate_Y_symmetric", (DL_FUNC) &_TestGLTM_calculate_Y_symmetric, 5},
    {"_TestGLTM_calculate_Y_m_dep", (DL_FUNC) &_TestGLTM_calculate_Y_m_dep, 4},
    {"_TestGLTM_bootstrap_independent", (DL_FUNC) &_TestGLTM_bootstrap_independent, 4},
    {"_TestGLTM_bootstrap_m_dep", (DL_FUNC) &_TestGLTM_bootstrap_m_dep, 6},
    {"_TestGLTM_bootstrap_U", (DL_FUNC) &_TestGLTM_bootstrap_U, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_TestGLTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
