#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double fourth_mom(NumericMatrix S, int a, int b, int c, int d) {
  double out;
  out = S(a,b) * S(c,d) + S(a,c) * S(b,d) + S(a,d) * S(b,c);
  return out;
}




// [[Rcpp::export]]
// NumericMatrix compute_cov_run_over(NumericMatrix S, IntegerMatrix ind_eq) {
//   
//   
//   // for (int j = 0; j < nr_ind_eq; j++) {
//   //   h_tilde[j] = X1[ind_eq(j,0)-1] *  X1[ind_eq(j,1)-1] * X2[ind_eq(j,2)-1] * X2[ind_eq(j,3)-1] 
//   //   - X1[ind_eq(j,4)-1] * X1[ind_eq(j,5)-1] * X2[ind_eq(j,6)-1] * X2[ind_eq(j,7)-1];
//   // }
//   
//   int nr_cols = ind_eq.nrow();
//   int p;
//   int q;
//   int r;
//   int s;
//   int u;
//   int v;
//   int w;
//   int z;
//   NumericMatrix cov(nr_cols,nr_cols);
//   
//   for (int i = 0; i < nr_cols; i++) {
//     p = indices(i,0)-1;
//     q = indices(i,1)-1;
//     r = indices(i,2)-1;
//     s = indices(i,3)-1;
//     for (int j = 0; j < nr_cols; j++) {
//       u = indices(j,0)-1;
//       v = indices(j,1)-1;
//       w = indices(j,2)-1;
//       z = indices(j,3)-1;
//       
//       cov(i,j) = 
//         fourth_mom(S, p, r, u, w) * fourth_mom(S, q, s, v, z) -
//         fourth_mom(S, p, r, u, z) * fourth_mom(S, q, s, v, w) -
//         fourth_mom(S, p, s, u, w) * fourth_mom(S, q, r, v, z) +
//         fourth_mom(S, p, s, u, z) * fourth_mom(S, q, r, v, w) +
//         S(p,r) * fourth_mom(S, q, s, u, w) * S(v,z) -
//         S(p,r) * fourth_mom(S, q, s, u, z) * S(v,w) -
//         S(p,s) * fourth_mom(S, q, r, u, w) * S(v,z) +
//         S(p,s) * fourth_mom(S, q, r, u, z) * S(v,w) +
//         S(u,w) * fourth_mom(S, p, r, v, z) * S(q,s) -
//         S(u,z) * fourth_mom(S, p, r, v, w) * S(q,s) -
//         S(u,w) * fourth_mom(S, p, s, v, z) * S(q,r) +
//         S(u,z) * fourth_mom(S, p, s, v, w) * S(q,r) -
//         3 * ((S(p,r) * S(q,s) - S(p,s) * S(q,r)) * (S(u,w) * S(v,z) - S(u,z) * S(v,w)));
//       
//       cov(j,i) = cov(i,j);
//     }
//   }
//   return cov;
// }
// 
