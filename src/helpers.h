#include <Rcpp.h>
#ifndef HELPERS
#define HELPERS

Rcpp::List random_combs(int n, int k, int nr);

int findn(int N, int D);

Rcpp::IntegerMatrix compute_S(int n, int i, int L);

Rcpp::IntegerMatrix permutations(int n);

Rcpp::List update_param(Rcpp::NumericMatrix S, 
                        Rcpp::NumericMatrix edges, 
                        int nr_obs);

Rcpp::NumericMatrix cov_from_graph_large(Rcpp::NumericVector Omega, 
                                         Rcpp::NumericVector Rho, 
                                         Rcpp::List paths);

int binomialCoefficient(int n, int k);

Rcpp::IntegerMatrix two_subsets(int n);

#endif