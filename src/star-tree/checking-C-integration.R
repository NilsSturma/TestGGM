library(MASS)
library(CombMSC)
library(microbenchmark)

library(Rcpp)
# https://teuder.github.io/rcpp4everyone_en/
# Important: Indices start with 0 in C++.
sourceCpp('rcpp-functions.cpp')



###############
# R functions #
###############

fourth_momR = function(S, ind){
  S[ind[1],ind[2]] * S[ind[3],ind[4]] + S[ind[1],ind[3]] * S[ind[2],ind[4]] + S[ind[1],ind[4]] * S[ind[2],ind[3]]
}

compute_cov_YR = function(S, indices){
  nr_cols = nrow(indices)
  # calculate covariance of Y
  cov = matrix(0, nrow=nr_cols, ncol=nr_cols)
  for (i in 1:nr_cols){
    p = indices[i,1]
    q = indices[i,2]
    r = indices[i,3]
    s = indices[i,4]
    for (j in i:nr_cols){
      u = indices[j,1]
      v = indices[j,2]
      w = indices[j,3]
      z = indices[j,4]

      cov[i,j] = 
        fourth_momR(S, c(p,r,u,w)) * fourth_momR(S, c(q,s,v,z)) - 
        fourth_momR(S, c(p,r,u,z)) * fourth_momR(S, c(q,s,v,w)) -
        fourth_momR(S, c(p,s,u,w)) * fourth_momR(S, c(q,r,v,z)) +
        fourth_momR(S, c(p,s,u,z)) * fourth_momR(S, c(q,r,v,w)) +
        (
          S[p,r] * fourth_momR(S, c(q,s,u,w)) * S[v,z] -
            S[p,r] * fourth_momR(S, c(q,s,u,z)) * S[v,w] -
            S[p,s] * fourth_momR(S, c(q,r,u,w)) * S[v,z] +
            S[p,s] * fourth_momR(S, c(q,r,u,z)) * S[v,w]
        ) +
        (
          S[u,w] * fourth_momR(S, c(p,r,v,z)) * S[q,s] -
            S[u,z] * fourth_momR(S, c(p,r,v,w)) * S[q,s] -
            S[u,w] * fourth_momR(S, c(p,s,v,z)) * S[q,r] +
            S[u,z] * fourth_momR(S, c(p,s,v,w)) * S[q,r]
        ) - 
        3 * (
          (S[p,r]*S[q,s] - S[p,s]*S[q,r]) * (S[u,w]*S[v,z] - S[u,z]*S[v,w])
        )
      cov[j,i] = cov[i,j]
    }
  }
  return(cov)
}

################
# benchmarking #
################

# generate data
n = 500
m = 10
beta = rep(1,m)
Sigma = beta %*% t(beta) + diag(rep(1,m))
X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma)

X1 = X[1:(n/2),]
X2 = X[((n/2)+1):n, ]

sub_sets = subsets(m,4,1:m)
nr_cols = 2*nrow(sub_sets)
indices = matrix(0, nrow=nr_cols, 4)
indices[1:(nr_cols/2),] = sub_sets
indices[((nr_cols/2)+1):(2*(nr_cols/2)),] = sub_sets[,c(1,3,2,4)]
mode(indices) = "integer"


# Sample covariance 
S = matrix(0, nrow=m, ncol=m)
for (i in 1:n){
  S = S + X[i,] %*% t(X[i,])
}
S = S/n

resR = compute_cov_YR(S, indices)
res = compute_cov_Y(S,indices)

microbenchmark(compute_cov_YR(S, indices),
               compute_cov_Y(S,indices),
               times=10L)
