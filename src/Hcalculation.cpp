#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector h_tilde(NumericVector X1,
                     NumericVector X2,
                     NumericVector X3,
                     NumericVector X4,
                     IntegerMatrix ind_eq, 
                     IntegerMatrix ind_ineq1, 
                     IntegerMatrix ind_ineq2){
  
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  NumericVector h_tilde(nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  ind_eq = ind_eq - 1;
  ind_ineq1 = ind_ineq1 - 1;
  ind_ineq2 = ind_ineq2 - 1;
  
  for (int j = 0; j < nr_ind_eq; j++) {
    h_tilde[j] = X1[ind_eq(j,0)] *  X1[ind_eq(j,1)] * X2[ind_eq(j,2)] * X2[ind_eq(j,3)] 
    - X1[ind_eq(j,4)] * X1[ind_eq(j,5)] * X2[ind_eq(j,6)] * X2[ind_eq(j,7)];
  }
  for (int j = 0; j < nr_ind_ineq1; j++) {
    h_tilde[nr_ind_eq + j] = - X1[ind_ineq1(j,0)] *  X1[ind_ineq1(j,1)] 
    * X2[ind_ineq1(j,2)] * X2[ind_ineq1(j,3)] 
    * X3[ind_ineq1(j,4)] * X3[ind_ineq1(j,5)];
  }
  for (int j = 0; j < nr_ind_ineq2; j++) {
    h_tilde[nr_ind_eq + nr_ind_ineq1 + j] =  X1[ind_ineq2(j,0)] * X1[ind_ineq2(j,1)] * X2[ind_ineq2(j,0)] * X2[ind_ineq2(j,1)] 
    * X3[ind_ineq2(j,2)] * X3[ind_ineq2(j,3)] * X4[ind_ineq2(j,2)] * X4[ind_ineq2(j,3)] 
    - X1[ind_ineq2(j,4)] * X1[ind_ineq2(j,5)] * X2[ind_ineq2(j,4)] * X2[ind_ineq2(j,5)] 
    * X3[ind_ineq2(j,6)] * X3[ind_ineq2(j,7)] * X4[ind_ineq2(j,6)] * X4[ind_ineq2(j,7)];
  }
  return(h_tilde);
}


// [[Rcpp::export]]
NumericVector h(List L, 
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2,
                IntegerMatrix perm){
  
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  NumericVector h(nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  
  for (int per = 0; per < perm.nrow(); per++){
    h = h + h_tilde(L[perm(per,0)], L[perm(per,1)], L[perm(per,2)], L[perm(per,3)], ind_eq, ind_ineq1, ind_ineq2);
  }
  return(h/perm.nrow());
}


// [[Rcpp::export]]
NumericMatrix calculate_H(NumericMatrix X,
                IntegerMatrix indices_U, 
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2){
  
  NumericMatrix H(indices_U.nrow(), ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  indices_U = indices_U - 1;
  
  IntegerMatrix perm = permutations(4);
  perm = perm - 1;

  for (int i = 0; i < indices_U.nrow(); i++){
    List L = List::create(X(indices_U(i,0),_), X(indices_U(i,1),_), X(indices_U(i,2),_), X(indices_U(i,3),_));
    H(i,_) = h(L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(H);
}


// [[Rcpp::export]]
NumericMatrix calculate_H_not_symmetric(NumericMatrix X,
                          IntegerMatrix indices_U, 
                          IntegerMatrix ind_eq, 
                          IntegerMatrix ind_ineq1, 
                          IntegerMatrix ind_ineq2){
  
  NumericMatrix H(indices_U.nrow(), ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  indices_U = indices_U - 1;
  for (int i = 0; i < indices_U.nrow(); i++){
    H(i,_) = h_tilde(X(indices_U(i,0),_), X(indices_U(i,1),_), X(indices_U(i,2),_), X(indices_U(i,3),_), ind_eq, ind_ineq1, ind_ineq2);
  }
  return(H);
}



// [[Rcpp::export]]
NumericVector calculate_g(NumericMatrix X,
                int i,
                int L,
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2,
                IntegerMatrix perm){
  
  int n = X.nrow();
  NumericVector g(ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  int K = floor((n-1)/L);
  
  IntegerMatrix S = compute_S(n,i,L);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    List L = List::create(X((i-1),_), X(S(k,0),_), X(S(k,1),_), X(S(k,2),_));
    g = g + h(L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(g/K);
}



// [[Rcpp::export]]
NumericMatrix calculate_G(NumericMatrix X,
                          int L,
                          IntegerMatrix ind_eq, 
                          IntegerMatrix ind_ineq1, 
                          IntegerMatrix ind_ineq2){
  
  int n = X.nrow();
  NumericMatrix G(n, ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  
  IntegerMatrix perm = permutations(4);
  perm = perm - 1;

  for (int i = 1; i <= n; i++){
    G((i-1),_) = calculate_g(X, i, L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(G);
}






/////////////////////////////////////////////////////
// calculate G faster with n_1 = K < sample_size n //
/////////////////////////////////////////////////////

// [[Rcpp::export]]
IntegerMatrix compute_S_small(int n, int i, int L, int K){
  
  IntegerVector v(n-1);
  for (int j = 1; j < (n+1); ++j){
    if (j < i){
      v[j-1] = j;
    } else if (j==i) {
      continue;
    } else {
      v[j-2] = j;
    }
  }
  // Only use n1 sets S
  int K_max = floor((n-1)/L);
  IntegerVector K_sets = sample(K_max,K)-1;
  IntegerVector sequence(L);
  IntegerMatrix res(K,L);
  
  for (int k = 0; k < K; ++k){
    sequence = seq((K_sets[k]*L), (((K_sets[k]+1)*L)-1));
    for (int l = 0; l < L; l++){
      res(k,l)=v[sequence[l]];
    }
  }
  return res;
}

// [[Rcpp::export]]
NumericVector g_small(NumericMatrix X,
                int i,
                int L,
                int n1,
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2,
                IntegerMatrix perm){
  
  int n = X.nrow();
  NumericVector g(ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  
  int K_max = floor((n-1)/L);
  int K =  min(IntegerVector::create(n1, K_max));
  
  IntegerMatrix S = compute_S_small(n,i,L,K);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    List L = List::create(X((i-1),_), X(S(k,0),_), X(S(k,1),_), X(S(k,2),_));
    g = g + h(L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(g/K);
}



// [[Rcpp::export]]
NumericMatrix calculate_G_small(NumericMatrix X,
                          int L,
                          int n1,
                          IntegerMatrix ind_eq, 
                          IntegerMatrix ind_ineq1, 
                          IntegerMatrix ind_ineq2){
  
  int n = X.nrow();
  NumericMatrix G(n1, ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  
  IntegerMatrix perm = permutations(4);
  perm = perm - 1;
  
  IntegerVector I = sample(n,n1);
  
  for (int i = 0; i < n1; i++){
    G(i,_) = g_small(X, I[i], L, n1, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(G);
}













/////////////////////////////
// Testing just equalities //
/////////////////////////////

// [[Rcpp::export]]
NumericVector h_tilde_eq(NumericVector X1,
                         NumericVector X2,
                         IntegerMatrix ind_eq){
  
  const int nr_ind_eq = ind_eq.nrow();
  NumericVector h_tilde(nr_ind_eq);
  ind_eq = ind_eq - 1;
  
  for (int j = 0; j < nr_ind_eq; j++) {
    h_tilde[j] = X1[ind_eq(j,0)] *  X1[ind_eq(j,1)] * X2[ind_eq(j,2)] * X2[ind_eq(j,3)] 
    - X1[ind_eq(j,4)] * X1[ind_eq(j,5)] * X2[ind_eq(j,6)] * X2[ind_eq(j,7)];
  }
  return(h_tilde);
}


// [[Rcpp::export]]
NumericVector h_eq(List L, 
                   IntegerMatrix ind_eq, 
                   IntegerMatrix perm){
  
  const int nr_ind_eq = ind_eq.nrow();
  NumericVector h(nr_ind_eq);
  
  for (int per = 0; per < perm.nrow(); per++){
    h = h + h_tilde_eq(L[perm(per,0)], L[perm(per,1)], ind_eq);
  }
  return(h/perm.nrow());
}


// [[Rcpp::export]]
NumericMatrix calculate_H_eq(NumericMatrix X,
                             IntegerMatrix indices_U, 
                             IntegerMatrix ind_eq){
  
  NumericMatrix H(indices_U.nrow(), ind_eq.nrow());
  indices_U = indices_U - 1;
  
  IntegerMatrix perm = permutations(2);
  perm = perm - 1;
  
  for (int i = 0; i < indices_U.nrow(); i++){
    List L = List::create(X(indices_U(i,0),_), X(indices_U(i,1),_));
    H(i,_) = h_eq(L, ind_eq, perm);
  }
  return(H);
}


// [[Rcpp::export]]
NumericMatrix calculate_H_not_symmetric_eq(NumericMatrix X,
                                           IntegerMatrix indices_U, 
                                           IntegerMatrix ind_eq){
  
  NumericMatrix H(indices_U.nrow(), ind_eq.nrow());
  indices_U = indices_U - 1;
  for (int i = 0; i < indices_U.nrow(); i++){
    H(i,_) = h_tilde_eq(X(indices_U(i,0),_), X(indices_U(i,1),_), ind_eq);
  }
  return(H);
}



// [[Rcpp::export]]
NumericVector g_eq(NumericMatrix X,
                   int i,
                   int L,
                   IntegerMatrix ind_eq, 
                   IntegerMatrix perm){
  
  int n = X.nrow();
  NumericVector g(ind_eq.nrow());
  int K = floor((n-1)/L);
  
  IntegerMatrix S = compute_S(n,i,L);
  S = S - 1;
  
  for (int k = 0; k < K; k++){
    List L = List::create(X((i-1),_), X(S(k,0),_));
    g = g + h_eq(L, ind_eq, perm);
  }
  return(g/K);
}


// [[Rcpp::export]]
NumericMatrix calculate_G_eq(NumericMatrix X,
                             int L,
                             IntegerMatrix ind_eq){
  
  int n = X.nrow();
  NumericMatrix G(n, ind_eq.nrow());
  
  IntegerMatrix perm = permutations(2);
  perm = perm - 1;
  
  for (int i = 1; i <= n; i++){
    G((i-1),_) = g_eq(X, i, L, ind_eq, perm);
  }
  return(G);
}



