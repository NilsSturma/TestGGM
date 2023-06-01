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
  
  while(sub_sets.size() < (unsigned)nr){
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




// [[Rcpp::export]]
List update_param(NumericMatrix S, NumericMatrix edges, int nr_obs){
  
  //update edge parameter
  NumericVector Rho(edges.nrow());
  for(int i=0; i < edges.nrow(); i++){
    Rho[i] = S((edges(i,0)-1),(edges(i,1)-1)) 
      / (sqrt(S((edges(i,0)-1),(edges(i,0)-1))) * sqrt(S((edges(i,1)-1),(edges(i,1)-1))));
  }
  
  //update observed variances
  NumericVector Omega(nr_obs);
  for(int i=0; i < nr_obs; i++){
    Omega[i] = S(i,i);
  }
  
  List res;
  res["Rho"] = Rho;
  res["Omega"] = Omega;
  
  return res;
}


// [[Rcpp::export]]
NumericMatrix cov_from_graph_large(NumericVector Omega, NumericVector Rho, List paths){
  
  // create vector of standard deviations
  int nr_obs = Omega.length();
  int nr_hidden = paths.length() - nr_obs;
  int nr_total = nr_obs+nr_hidden;
  NumericVector std(nr_total);
  for (int i=0; i < nr_obs; i++){
    std[i] = sqrt(Omega[i]);
  }
  for (int i=0; i < nr_hidden; i++){
    std[nr_obs+i] = 1;
  }
  
  // create large covariance matrix
  double prod = 1;
  NumericMatrix cov(nr_total, nr_total);
  
  for (int i=0; i < nr_total; i++){
    cov(i,i) = std[i] * std[i];
  }
  for (int i=0; i < (nr_total-1); i++){
    List paths_from_i = paths[i];
    for(int j=1; j < nr_total; j++){
      prod = 1;
      IntegerVector path = paths_from_i[j];  //this vector has always a different length
      for (int k=0; k < path.length(); k++){
        prod = prod * Rho[(path[k]-1)];
      }
      cov(i,j) = std[i] * std[j] * prod;
      cov(j,i) = cov(i,j);
    }
  }
  return cov;
}

// [[Rcpp::export]]
int binomialCoefficient(int n, int k) {
  if (k > n - k)
    k = n - k;
  int res = 1;
  for (int i = 0; i < k; ++i) {
    res *= (n - i);
    res /= (i + 1);
  }
  return res;
}

// [[Rcpp::export]]
IntegerMatrix two_subsets(int n){
  IntegerMatrix res(binomialCoefficient(n,2),2);
  int count = 0;
  for(int i=1; i <= n; i++){
    for(int j=i+1; j <= n; j++){
      res(count,0) = i;
      res(count,1) = j;
      count = count + 1;
    }
  }
  return res;
}

