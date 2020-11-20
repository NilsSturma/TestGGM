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
  
  const int n = findn(X.nrow(),4)/4;
  const NumericMatrix X1 = X(Range(0,(n-1)),_);
  const NumericMatrix X2 = X(Range(n,(2*n-1)),_);
  const NumericMatrix X3 = X(Range(2*n,(3*n-1)),_);
  const NumericMatrix X4 = X(Range(3*n,(4*n-1)),_);
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  
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
NumericMatrix calculate_Y_symmetric(NumericMatrix X, 
                                    IntegerMatrix ind_eq, 
                                    IntegerMatrix ind_ineq1, 
                                    IntegerMatrix ind_ineq2, 
                                    IntegerMatrix perm) {
  
  perm = perm - 1;
  const int n = findn(X.nrow(),4)/4;
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  NumericMatrix Y(n, nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  
  NumericMatrix X1(n,X.ncol()); 
  NumericMatrix X2(n,X.ncol()); 
  NumericMatrix X3(n,X.ncol()); 
  NumericMatrix X4(n,X.ncol()); 
  
  
  for (int per = 0; per < perm.nrow(); per++){
    
    X1 = X(Range(perm(per,0)*n,((perm(per,0)+1)*n-1)),_);
    X2 = X(Range(perm(per,1)*n,((perm(per,1)+1)*n-1)),_);
    X3 = X(Range(perm(per,2)*n,((perm(per,2)+1)*n-1)),_);
    X4 = X(Range(perm(per,3)*n,((perm(per,3)+1)*n-1)),_);
  
    for (int j = 0; j < nr_ind_eq; j++) {
      Y(_,j) = Y(_,j) + X1(_,ind_eq(j,0)-1) *  X1(_,ind_eq(j,1)-1) * X2(_,ind_eq(j,2)-1) * X2(_,ind_eq(j,3)-1) 
      - X1(_,ind_eq(j,4)-1) * X1(_,ind_eq(j,5)-1) * X2(_,ind_eq(j,6)-1) * X2(_,ind_eq(j,7)-1);
    }
    
    for (int j = 0; j < nr_ind_ineq1; j++) {
      Y(_,nr_ind_eq + j) = Y(_,nr_ind_eq + j) - X1(_,ind_ineq1(j,0)-1) *  X1(_,ind_ineq1(j,1)-1) 
      * X2(_,ind_ineq1(j,2)-1) * X2(_,ind_ineq1(j,3)-1) 
      * X3(_,ind_ineq1(j,4)-1) * X3(_,ind_ineq1(j,5)-1);
    }
    
    for (int j = 0; j < nr_ind_ineq2; j++) {
      Y(_,nr_ind_eq + nr_ind_ineq1 + j) = Y(_,nr_ind_eq + nr_ind_ineq1 + j) + X1(_,ind_ineq2(j,0)-1) * X1(_,ind_ineq2(j,1)-1) * X2(_,ind_ineq2(j,0)-1) * X2(_,ind_ineq2(j,1)-1) 
      * X3(_,ind_ineq2(j,2)-1) * X3(_,ind_ineq2(j,3)-1) * X4(_,ind_ineq2(j,2)-1) * X4(_,ind_ineq2(j,3)-1) 
      - X1(_,ind_ineq2(j,4)-1) * X1(_,ind_ineq2(j,5)-1) * X2(_,ind_ineq2(j,4)-1) * X2(_,ind_ineq2(j,5)-1) 
      * X3(_,ind_ineq2(j,6)-1) * X3(_,ind_ineq2(j,7)-1) * X4(_,ind_ineq2(j,6)-1) * X4(_,ind_ineq2(j,7)-1);
    }
  } 
  return (Y/perm.nrow());
  
}







// [[Rcpp::export]]
NumericMatrix calculate_Y_m_dep(NumericMatrix X, IntegerMatrix ind_eq, IntegerMatrix ind_ineq1, IntegerMatrix ind_ineq2) {
  
  const int n = X.nrow();
  const NumericMatrix X1 = X(Range(0,(n-4)),_);
  const NumericMatrix X2 = X(Range(1,(n-3)),_);
  const NumericMatrix X3 = X(Range(2,(n-2)),_);
  const NumericMatrix X4 = X(Range(3,(n-1)),_);
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  
  NumericMatrix Y((n-3), nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  
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