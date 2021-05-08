// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <set>
using namespace Rcpp;


// [[Rcpp::export]]
List random_combs(int n, int k, int nr){
  std::set<std::set<int>> sub_sets;
  
  IntegerVector v(n);
  
  for(int i=0; i < n; i++){
    v[i] = i+1;
  }
  
  std::set<int> s;
  IntegerVector w(k);
  
  while(sub_sets.size() < nr){
    w = sample(v,k);
    for (int i = 0; i < k; i++) {
      s.insert(w[i]);
    }
    sub_sets.insert(s);
    s.clear();
  }
  
  List res = List::create(sub_sets);
  
  return(res);
}



// [[Rcpp::export]]
int findn(int N, int D){ 
  int rem = N % D; 
  
  if (rem == 0) 
    return N; 
  else
    return N - rem; 
} 



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

