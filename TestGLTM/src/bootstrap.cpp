#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix bootstrap_independent(int E, NumericMatrix Y_centered){
  
  const int n = Y_centered.nrow();
  const int p = Y_centered.ncol();
  
  NumericVector epsilons(n);
  NumericVector colsums(p);
  NumericMatrix W(E,p);
  
  for (int i = 0; i < E; i++){
    epsilons = rnorm(n,0,1);
    for (int j = 0; j < p; j++){
      colsums[j] = sum(Y_centered(_,j) * epsilons);
    }
    W(i,_) = (1/sqrt(n)) * colsums;
  }
  return W;
}





// [[Rcpp::export]]
NumericMatrix bootstrap_m_dep(int E, int B, int omega, NumericMatrix Y_centered){
  
  const int p = Y_centered.ncol();
  
  NumericVector epsilons(omega);
  double batchsum = 0;
  NumericVector sums(p);
  IntegerVector L(B);
  NumericMatrix W(E,p);
  
  for (int i = 0; i < E; i++){
    epsilons = rnorm(omega,0,1);
    sums.fill(0);
    for (int j = 0; j < p; j++){
      for (int b = 0; b < omega; b++){
        L = seq((b*B), (((b+1)*B)-1));
        batchsum = 0;
        for (int l = 0; l < B; l++){
          batchsum = batchsum + Y_centered(L[l],j);
        }
        sums[j] = sums[j] + batchsum * epsilons[b];
      }
    }
    W(i,_) = (1/sqrt(B * omega)) * sums;
  }
  return W;
}


// [[Rcpp::export]]
List bootstrap_U(int E, NumericMatrix H_centered, NumericMatrix G_centered, int N){
  
  int N_hat = H_centered.nrow();
  int n = G_centered.nrow();
  int p = H_centered.ncol();
  
  NumericVector epsilons((N_hat + n));
  NumericVector colsums_B(p);
  NumericMatrix U_B(E,p);
  NumericVector colsums_A(p);
  NumericMatrix U_A(E,p);
  
  for (int i = 0; i < E; i++){
    epsilons = rnorm((N_hat+n),0,1);
    for (int j = 0; j < p; j++){
      colsums_B[j] = sum(H_centered(_,j) * epsilons[Range(0, (N_hat-1))]);
    }
    U_B(i,_) = (1/sqrt(N_hat)) * colsums_B;
    for (int j = 0; j < p; j++){
      colsums_A[j] = sum(G_centered(_,j) * epsilons[Range(N_hat, ((N_hat+n)-1))]);
    }
    U_A(i,_) = (4/sqrt(n)) * colsums_A;
  }
  List res = List::create(U_A, U_B);
  return res;
}