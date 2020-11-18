#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
int findn(int N, int D){ 
  int rem = N % D; 
  
  if (rem == 0) 
    return N; 
  else
    return N - rem; 
} 



// [[Rcpp::export]]
NumericMatrix calculate_Y_independent(NumericMatrix X, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2) {
  
  int n;
  n = findn(X.nrow(),4)/4;
  
  NumericMatrix X1 = X(Range(0,(n-1)),_);
  NumericMatrix X2 = X(Range(n,(2*n-1)),_);
  NumericMatrix X3 = X(Range(2*n,(3*n-1)),_);
  NumericMatrix X4 = X(Range(3*n,(4*n-1)),_);
  
  int nr_ind_eq = ind_eq.nrow();
  int nr_ind_ineq1 = ind_ineq1.nrow();
  int nr_ind_ineq2 = ind_ineq2.nrow();
  
  NumericMatrix Y(n, nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);

  for (int j = 0; j < nr_ind_eq; j++) {
    Y(_,j) = X1(_,ind_eq(j,0)-1) *  X1(_,ind_eq(j,1)-1) * X2(_,ind_eq(j,2)-1) * X2(_,ind_eq(j,3)-1) 
      - X1(_,ind_eq(j,4)-1) * X1(_,ind_eq(j,5)-1) * X2(_,ind_eq(j,6)-1) * X2(_,ind_eq(j,7)-1);
  }
  
  for (int j = 0; j < nr_ind_ineq1; j++) {
    Y(_,nr_ind_eq + j) = - X1(_,ind_ineq1(j,0)-1) *  X1(_,ind_ineq1(j,1)-1) 
      * X2(_,ind_ineq1(j,2)-1) * X2(_,ind_ineq1(j,3)-1) 
      * X3(_,ind_ineq1(j,4)-1) * X3(_,ind_ineq1(j,5)-1);
  }
  
  for (int j = 0; j < nr_ind_ineq2; j++) {
    Y(_,nr_ind_eq + nr_ind_ineq1 + j) = X1(_,ind_ineq2(j,0)-1) * X1(_,ind_ineq2(j,1)-1) * X2(_,ind_ineq2(j,0)-1) * X2(_,ind_ineq2(j,1)-1) 
    * X3(_,ind_ineq2(j,2)-1) * X3(_,ind_ineq2(j,3)-1) * X4(_,ind_ineq2(j,2)-1) * X4(_,ind_ineq2(j,3)-1) 
    - X1(_,ind_ineq2(j,4)-1) * X1(_,ind_ineq2(j,5)-1) * X2(_,ind_ineq2(j,4)-1) * X2(_,ind_ineq2(j,5)-1) 
    * X3(_,ind_ineq2(j,6)-1) * X3(_,ind_ineq2(j,7)-1) * X4(_,ind_ineq2(j,6)-1) * X4(_,ind_ineq2(j,7)-1);
  }
  
  return Y;
  
}



// [[Rcpp::export]]
NumericVector bootstrap_independent(int E, NumericVector standardizer, NumericMatrix Y_centered, int p_eq){
  
  int n = Y_centered.nrow();
  int p = Y_centered.ncol();
  NumericVector epsilons(n);
  NumericVector colsums(p);
  NumericVector m(p);
  NumericVector res(E);
  
  for (int i = 0; i < E; i++){
    epsilons = rnorm(n,0,1);
    for (int j = 0; j < p; j++){
      colsums[j] = sum(Y_centered(_,j) * epsilons);
    }
    m[Range(0, (p_eq-1))] = abs(standardizer[Range(0, (p_eq-1))] * colsums[Range(0, (p_eq-1))]);
    m[Range(p_eq, (p-1))] = standardizer[Range(p_eq, (p-1))] * colsums[Range(p_eq, (p-1))];
    res[i] = (1/sqrt(n)) * max(m);
  }

  return res;
}
