#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector h_tilde_fac(NumericVector X1,
                          NumericVector X2,
                          NumericVector X3,
                          IntegerMatrix ind){
  ind = ind - 1;
  NumericVector h_tilde(ind.nrow());
  
  for (int j = 0; j < ind.nrow(); j++) {
    h_tilde[j] = X1[ind(j,0)] *  X1[ind(j,3)] * X2[ind(j,1)] * X2[ind(j,4)] * X3[ind(j,2)] * X3[ind(j,5)]
               + X1[ind(j,0)] *  X1[ind(j,4)] * X2[ind(j,1)] * X2[ind(j,5)] * X3[ind(j,2)] * X3[ind(j,3)]
               + X1[ind(j,0)] *  X1[ind(j,5)] * X2[ind(j,1)] * X2[ind(j,3)] * X3[ind(j,2)] * X3[ind(j,4)] 
               - X3[ind(j,2)] *  X3[ind(j,3)] * X2[ind(j,1)] * X2[ind(j,4)] * X1[ind(j,0)] * X1[ind(j,5)] 
               - X3[ind(j,2)] *  X3[ind(j,4)] * X2[ind(j,1)] * X2[ind(j,5)] * X1[ind(j,0)] * X1[ind(j,3)] 
               - X3[ind(j,2)] *  X3[ind(j,5)] * X2[ind(j,1)] * X2[ind(j,3)] * X1[ind(j,0)] * X1[ind(j,4)];
  }
  return(h_tilde);
}


// [[Rcpp::export]]
NumericVector h_fac(List L, 
                    IntegerMatrix ind_minors,
                    IntegerMatrix perm){
  
  NumericVector h(ind_minors.nrow());
  
  for (int per = 0; per < perm.nrow(); per++){
    h = h + h_tilde_fac(L[perm(per,0)], L[perm(per,1)], L[perm(per,2)], ind_minors);
  }
  return(h/perm.nrow());
}


// [[Rcpp::export]]
NumericMatrix H_factors(NumericMatrix X,
                        IntegerMatrix indices, 
                        IntegerMatrix ind_minors){
  
  NumericMatrix H(indices.nrow(), ind_minors.nrow());
  indices = indices - 1;
  
  IntegerMatrix perm = permutations((0.5*ind_minors.ncol()));
  perm = perm - 1;
  
  for (int i = 0; i < indices.nrow(); i++){
    List L = List::create(X(indices(i,0),_), X(indices(i,1),_), X(indices(i,2),_));
    H(i,_) = h_fac(L, ind_minors, perm);
  }
  return(H);
}



// [[Rcpp::export]]
NumericVector g_fac(NumericMatrix X,
                    int i,
                    int L,
                    IntegerMatrix ind,
                    IntegerMatrix perm){
  
  int n = X.nrow();
  NumericVector g(ind.nrow());
  int K = floor((n-1)/L);
  
  IntegerMatrix S = compute_S(n,i,L);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    List L = List::create(X((i-1),_), X(S(k,0),_), X(S(k,1),_));
    g = g + h_fac(L, ind, perm);
  }
  return(g/K);
}



// [[Rcpp::export]]
NumericMatrix G_factors(NumericMatrix X,
                          int L,
                          IntegerMatrix ind_minors){
  
  int n = X.nrow();
  NumericMatrix G(n, ind_minors.nrow());
  
  IntegerMatrix perm = permutations((0.5*ind_minors.ncol()));
  perm = perm - 1;
  
  for (int i = 1; i <= n; i++){
    G((i-1),_) = g_fac(X, i, L, ind_minors, perm);
  }
  return(G);
}
