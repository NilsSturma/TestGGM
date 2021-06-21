#include <Rcpp.h>
#ifndef HELPERS
#define HELPERS

Rcpp::List random_combs(int n, int k, int nr);

int findn(int N, int D);

Rcpp::IntegerMatrix compute_S(int n, int i, int L);

Rcpp::IntegerMatrix permutations(int n);

#endif