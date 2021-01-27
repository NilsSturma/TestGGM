#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double expectation(NumericMatrix S, IntegerVector j){
  double out;
  out = S(j[0], j[1]) * S(j[2], j[3]) - S(j[4], j[5]) * S(j[6], j[7]);
  return out;
}


// [[Rcpp::export]]
double fourth_mom(NumericMatrix S, int a, int b, int c, int d) {
  double out;
  out = S(a,b) * S(c,d) + S(a,c) * S(b,d) + S(a,d) * S(b,c);
  return out;
}

// [[Rcpp::export]]
NumericMatrix cov_grouping(NumericMatrix S, IntegerMatrix ind_eq){
  
  ind_eq = ind_eq - 1;
  NumericMatrix cov(ind_eq.nrow(),ind_eq.nrow());
  
  for (int j = 0; j < ind_eq.nrow(); j++){
    for (int k = j; k < ind_eq.nrow(); k++){
      cov(j,k) = 
        fourth_mom(S, ind_eq(j,0), ind_eq(j,1), ind_eq(k,0), ind_eq(k,1)) 
          * fourth_mom(S, ind_eq(j,2), ind_eq(j,3), ind_eq(k,2), ind_eq(k,3)) -
        fourth_mom(S, ind_eq(j,0), ind_eq(j,1), ind_eq(k,4), ind_eq(k,5)) 
          * fourth_mom(S, ind_eq(j,2), ind_eq(j,3), ind_eq(k,6), ind_eq(k,7)) -
        fourth_mom(S, ind_eq(j,4), ind_eq(j,5), ind_eq(k,0), ind_eq(k,1)) 
          * fourth_mom(S, ind_eq(j,6), ind_eq(j,7), ind_eq(k,2), ind_eq(k,3)) +
        fourth_mom(S, ind_eq(j,4), ind_eq(j,5), ind_eq(k,4), ind_eq(k,5)) 
          * fourth_mom(S, ind_eq(j,6), ind_eq(j,7), ind_eq(k,6), ind_eq(k,7)) -
        expectation(S,ind_eq(j,_)) * expectation(S,ind_eq(k,_));
      if (k != j){
        cov(k,j) = cov(j,k);
      }
    }
  }
  return cov;
}

// [[Rcpp::export]]
NumericMatrix cov_run_over(NumericMatrix S, IntegerMatrix ind_eq){
  
  NumericMatrix cov(ind_eq.nrow(),ind_eq.nrow());
  cov = cov_grouping(S, ind_eq);
  ind_eq = ind_eq - 1;
  
  for (int j = 0; j < ind_eq.nrow(); j++){
    for (int k = j; k < ind_eq.nrow(); k++){
      cov(j,k) = cov(j,k) +
        S(ind_eq(j,0), ind_eq(j,1)) * S(ind_eq(k,2), ind_eq(k,3)) * fourth_mom(S, ind_eq(j,2), ind_eq(j,3), ind_eq(k,0), ind_eq(k,1)) -
        S(ind_eq(j,0), ind_eq(j,1)) * S(ind_eq(k,6), ind_eq(k,7)) * fourth_mom(S, ind_eq(j,2), ind_eq(j,3), ind_eq(k,4), ind_eq(k,5)) -
        S(ind_eq(j,4), ind_eq(j,5)) * S(ind_eq(k,2), ind_eq(k,3)) * fourth_mom(S, ind_eq(j,6), ind_eq(j,7), ind_eq(k,0), ind_eq(k,1)) +
        S(ind_eq(j,4), ind_eq(j,5)) * S(ind_eq(k,6), ind_eq(k,7)) * fourth_mom(S, ind_eq(j,6), ind_eq(j,7), ind_eq(k,4), ind_eq(k,5)) + 
        S(ind_eq(k,0), ind_eq(k,1)) * S(ind_eq(j,2), ind_eq(j,3)) * fourth_mom(S, ind_eq(k,2), ind_eq(k,3), ind_eq(j,0), ind_eq(j,1)) -
        S(ind_eq(k,0), ind_eq(k,1)) * S(ind_eq(j,6), ind_eq(j,7)) * fourth_mom(S, ind_eq(k,2), ind_eq(k,3), ind_eq(j,4), ind_eq(j,5)) -
        S(ind_eq(k,4), ind_eq(k,5)) * S(ind_eq(j,2), ind_eq(j,3)) * fourth_mom(S, ind_eq(k,6), ind_eq(k,7), ind_eq(j,0), ind_eq(j,1)) +
        S(ind_eq(k,4), ind_eq(k,5)) * S(ind_eq(j,6), ind_eq(j,7)) * fourth_mom(S, ind_eq(k,6), ind_eq(k,7), ind_eq(j,4), ind_eq(j,5)) -
        2 * expectation(S,ind_eq(j,_)) * expectation(S,ind_eq(k,_));
      if (k != j){
        cov(k,j) = cov(j,k);
      }
    }
  }
  return cov;
}
