#include <Rcpp.h>
using namespace Rcpp;


//////////////////////
// helper functions //
//////////////////////

// [[Rcpp::export]]
IntegerMatrix compute_S(int n, int i, int L){
  int K = floor((n-1)/L);
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
  IntegerVector sequence(L);
  IntegerMatrix res(K,L);
  for (int k = 0; k < K; ++k){
    sequence = seq((k*L), (((k+1)*L)-1));
    for (int l = 0; l < L; l++){
      res(k,l)=v[sequence[l]];
    }
  }
  return res;
}


// [[Rcpp::export]]
IntegerMatrix permutations(int n){
  
  int factorial = 1;
  for (int i=1; i <= n; i++){
    factorial = factorial * i;
  }
  
  IntegerVector v(n);
  for (int i=1; i <=n; i++){
    v[i-1] = i;
  }
  
  IntegerMatrix res(factorial,n);
  for (int i=0; i<factorial; i++){
    res(i,_) = v; 
    std::next_permutation(v.begin(), v.end());
  }
  return res;
}


////////////////////
// main functions //
////////////////////

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
  
  for (int j = 0; j < nr_ind_eq; j++) {
    h_tilde[j] = X1[ind_eq(j,0)-1] *  X1[ind_eq(j,1)-1] * X2[ind_eq(j,2)-1] * X2[ind_eq(j,3)-1] 
    - X1[ind_eq(j,4)-1] * X1[ind_eq(j,5)-1] * X2[ind_eq(j,6)-1] * X2[ind_eq(j,7)-1];
  }
  for (int j = 0; j < nr_ind_ineq1; j++) {
    h_tilde[nr_ind_eq + j] = - X1[ind_ineq1(j,0)-1] *  X1[ind_ineq1(j,1)-1] 
    * X2[ind_ineq1(j,2)-1] * X2[ind_ineq1(j,3)-1] 
    * X3[ind_ineq1(j,4)-1] * X3[ind_ineq1(j,5)-1];
  }
  for (int j = 0; j < nr_ind_ineq2; j++) {
    h_tilde[nr_ind_eq + nr_ind_ineq1 + j] =  X1[ind_ineq2(j,0)-1] * X1[ind_ineq2(j,1)-1] * X2[ind_ineq2(j,0)-1] * X2[ind_ineq2(j,1)-1] 
    * X3[ind_ineq2(j,2)-1] * X3[ind_ineq2(j,3)-1] * X4[ind_ineq2(j,2)-1] * X4[ind_ineq2(j,3)-1] 
    - X1[ind_ineq2(j,4)-1] * X1[ind_ineq2(j,5)-1] * X2[ind_ineq2(j,4)-1] * X2[ind_ineq2(j,5)-1] 
    * X3[ind_ineq2(j,6)-1] * X3[ind_ineq2(j,7)-1] * X4[ind_ineq2(j,6)-1] * X4[ind_ineq2(j,7)-1];
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


// Calculation of g_i(X_i)
// [[Rcpp::export]]
NumericVector g(NumericMatrix X,
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
    List L = List::create(X(i,_), X(S(k,0),_), X(S(k,1),_), X(S(k,2),_));
    g = g + h(L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(g/K);
}


// Save all g_i(X_i) in matrix G
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

  for (int i = 0; i < n; i++){
    G(i,_) = g(X, i, L, ind_eq, ind_ineq1, ind_ineq2, perm);
  }
  return(G);
}
