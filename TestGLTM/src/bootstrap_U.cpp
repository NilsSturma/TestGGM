#include <Rcpp.h>
using namespace Rcpp;

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
