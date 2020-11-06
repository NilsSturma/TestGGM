#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix calculate_Y(IntegerMatrix indices, NumericMatrix X1, NumericMatrix X2) {
  int nr_cols = indices.nrow();
  int n = X1.nrow();
  int p = 0;
  int q = 0;
  int r = 0;
  int s = 0;
  NumericMatrix Y(n,nr_cols);
  for (int j = 0; j < nr_cols; j++) {
    p = indices(j,0)-1;
    q = indices(j,1)-1;
    r = indices(j,2)-1;
    s = indices(j,3)-1;
    Y(_,j) = X1(_,p) * X1(_,s) * X2(_,q) * X2(_,r) - X1(_,p) * X1(_,r) * X2(_,q) * X2(_,s);
  }
  return Y;
}


// [[Rcpp::export]]
double fourth_mom(NumericMatrix S, int a, int b, int c, int d) {
  double out;
  out = S(a,b) * S(c,d) + S(a,c) * S(b,d) + S(a,d) * S(b,c);
  return out;
}


// [[Rcpp::export]]
NumericMatrix compute_cov_Y(NumericMatrix S, IntegerMatrix indices) {
  int nr_cols = indices.nrow();
  int p;
  int q;
  int r;
  int s;
  int u;
  int v;
  int w;
  int z;
  NumericMatrix cov(nr_cols,nr_cols);
  
  for (int i = 0; i < nr_cols; i++) {
    p = indices(i,0)-1;
    q = indices(i,1)-1;
    r = indices(i,2)-1;
    s = indices(i,3)-1;
    for (int j = 0; j < nr_cols; j++) {
      u = indices(j,0)-1;
      v = indices(j,1)-1;
      w = indices(j,2)-1;
      z = indices(j,3)-1;
      
      cov(i,j) = 
        fourth_mom(S, p, r, u, w) * fourth_mom(S, q, s, v, z) -
        fourth_mom(S, p, r, u, z) * fourth_mom(S, q, s, v, w) -
        fourth_mom(S, p, s, u, w) * fourth_mom(S, q, r, v, z) +
        fourth_mom(S, p, s, u, z) * fourth_mom(S, q, r, v, w) +
        S(p,r) * fourth_mom(S, q, s, u, w) * S(v,z) -
        S(p,r) * fourth_mom(S, q, s, u, z) * S(v,w) -
        S(p,s) * fourth_mom(S, q, r, u, w) * S(v,z) +
        S(p,s) * fourth_mom(S, q, r, u, z) * S(v,w) +
        S(u,w) * fourth_mom(S, p, r, v, z) * S(q,s) -
        S(u,z) * fourth_mom(S, p, r, v, w) * S(q,s) -
        S(u,w) * fourth_mom(S, p, s, v, z) * S(q,r) +
        S(u,z) * fourth_mom(S, p, s, v, w) * S(q,r) -
        3 * ((S(p,r) * S(q,s) - S(p,s) * S(q,r)) * (S(u,w) * S(v,z) - S(u,z) * S(v,w)));
      
      cov(j,i) = cov(i,j);
    }
  }
  return cov;
}