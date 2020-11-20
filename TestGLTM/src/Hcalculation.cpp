#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector h_tile(NumericVector X1,
                     NumericVector X2,
                     NumericVector X3,
                     NumericVector X4,
                     IntegerMatrix ind_eq, 
                     IntegerMatrix ind_ineq1, 
                     IntegerMatrix ind_ineq2){
  
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  NumericVector h_tile(nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  
  for (int j = 0; j < nr_ind_eq; j++) {
    h_tile[j] = X1[ind_eq(j,0)-1] *  X1[ind_eq(j,1)-1] * X2[ind_eq(j,2)-1] * X2[ind_eq(j,3)-1] 
    - X1[ind_eq(j,4)-1] * X1[ind_eq(j,5)-1] * X2[ind_eq(j,6)-1] * X2[ind_eq(j,7)-1];
  }
  for (int j = 0; j < nr_ind_ineq1; j++) {
    h_tile[nr_ind_eq + j] = - X1[ind_ineq1(j,0)-1] *  X1[ind_ineq1(j,1)-1] 
    * X2[ind_ineq1(j,2)-1] * X2[ind_ineq1(j,3)-1] 
    * X3[ind_ineq1(j,4)-1] * X3[ind_ineq1(j,5)-1];
  }
  for (int j = 0; j < nr_ind_ineq2; j++) {
    h_tile[nr_ind_eq + nr_ind_ineq1 + j] =  X1[ind_ineq2(j,0)-1] * X1[ind_ineq2(j,1)-1] * X2[ind_ineq2(j,0)-1] * X2[ind_ineq2(j,1)-1] 
    * X3[ind_ineq2(j,2)-1] * X3[ind_ineq2(j,3)-1] * X4[ind_ineq2(j,2)-1] * X4[ind_ineq2(j,3)-1] 
    - X1[ind_ineq2(j,4)-1] * X1[ind_ineq2(j,5)-1] * X2[ind_ineq2(j,4)-1] * X2[ind_ineq2(j,5)-1] 
    * X3[ind_ineq2(j,6)-1] * X3[ind_ineq2(j,7)-1] * X4[ind_ineq2(j,6)-1] * X4[ind_ineq2(j,7)-1];
  }
  return(h_tile);
}



// [[Rcpp::export]]
NumericVector h(List L, 
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2){
  
  const int nr_ind_eq = ind_eq.nrow();
  const int nr_ind_ineq1 = ind_ineq1.nrow();
  const int nr_ind_ineq2 = ind_ineq2.nrow();
  NumericVector h(nr_ind_eq+nr_ind_ineq1+nr_ind_ineq2);
  
  Environment myEnv = Environment::global_env();
  Function perms = myEnv["permutations"];
  IntegerMatrix perm = as<IntegerMatrix>(perms(Rcpp::Named("n", L.length())));
  
  perm = perm - 1;
  
  for (int per = 0; per < perm.nrow(); per++){
    h = h + h_tile(L[perm(per,0)], L[perm(per,1)], L[perm(per,2)], L[perm(per,3)], ind_eq, ind_ineq1, ind_ineq2);
  }
  return(h/perm.nrow());
}


NumericMatrix H(NumericMatrix X,
                IntegerMatrix indices_U, 
                IntegerMatrix ind_eq, 
                IntegerMatrix ind_ineq1, 
                IntegerMatrix ind_ineq2){
  
  NumericMatrix H(indices_U.nrow(), ind_eq.nrow()+ind_ineq1.nrow()+ind_ineq2.nrow());
  List L = List::create(X(0,_), X(1,_), X(2,_), X(3,_));
  
  indices_U = indices_U - 1;
  
  for (int i = 0; i < indices_U.nrow(); i++){
    L = List::create(X(indices_U(i,0),_), X(indices_U(i,1),_), X(indices_U(i,2),_), X(indices_U(i,3),_));
    H(i,_) = h(L, ind_eq, ind_ineq1, ind_ineq2);
  }
  return(H);
}