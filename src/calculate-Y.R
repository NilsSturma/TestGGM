library(MASS) #mvrnorm
library(CombMSC) #subsets



# Function to compute fourth moments (covariance matrix given)
fourth_mom = function(S, ind){
  # S: Covariance matrix
  # ind: Vector of length 4. Indices to compute the fourht moment.  Example: ind=c(1,2,3,4) -> compute moment E[X_1*X_2*X_3*X_4]
  S[ind[1],ind[2]] * S[ind[3],ind[4]] + S[ind[1],ind[3]] * S[ind[2],ind[4]] + S[ind[1],ind[4]] * S[ind[2],ind[3]]
}




# Function to calculate Y ((n-1)x(2*choose(m,4))-matrix)
# Just for star tree!
compute_Y = function(X){
  m = dim(X)[2]
  n = dim(X)[1]
  
  # create all indices in right order !!!different order than previously!!!
  sub_sets = subsets(m,4,1:m)
  indices = matrix(0, nrow=2*nrow(sub_sets), 4)
  indices[1:nrow(sub_sets),] = sub_sets
  indices[(nrow(sub_sets)+1):(2*nrow(sub_sets)),] = sub_sets[,c(1,3,2,4)]
  
  # Compute Y (unbiased estimate of the tetrads)
  Y = matrix(0, nrow = n-1, ncol=(2*choose(m,4)))
  
  for (i in 1:(2*choose(m,4))){
    p = indices[i,1]
    q = indices[i,2]
    r = indices[i,3]
    s = indices[i,4]
    Y[,i] = X[1:(n-1),p] * X[1:(n-1),r] * X[2:n,q] * X[2:n,s] - X[1:(n-1),p] * X[1:(n-1),s] * X[2:n,q] * X[2:n,r]
  }
  return(Y)
}



# Function to calculate covariance of Y, usind the theoretical covariance matrix sigma
# Just for star tree!
compute_cov_Y = function(X, S){
  
  m = dim(X)[2]
  n = dim(X)[1]
  
  
  # create all indices in right order !!!different order than previously!!!
  sub_sets = subsets(m,4,1:m)
  indices = matrix(0, nrow=2*nrow(sub_sets), 4)
  indices[1:nrow(sub_sets),] = sub_sets
  indices[(nrow(sub_sets)+1):(2*nrow(sub_sets)),] = sub_sets[,c(1,3,2,4)]
  
  # calculate covariance of Y
  cov = matrix(0, nrow=(2*choose(m,4)), ncol=(2*choose(m,4)))
  # !!! could exploit symmetry here !!!
  for (i in 1:(2*choose(m,4))){
    p = indices[i,1]
    q = indices[i,2]
    r = indices[i,3]
    s = indices[i,4]
    for (j in 1:(2*choose(m,4))){
      u = indices[j,1]
      v = indices[j,2]
      w = indices[j,3]
      z = indices[j,4]
      cov[i,j] = 
        fourth_mom(S, c(p,r,u,w)) * fourth_mom(S, c(q,s,v,z)) - 
        fourth_mom(S, c(p,r,u,z)) * fourth_mom(S, c(q,s,v,w)) -
        fourth_mom(S, c(p,s,u,w)) * fourth_mom(S, c(q,r,v,z)) +
        fourth_mom(S, c(p,s,u,z)) * fourth_mom(S, c(q,r,v,w)) +
        (
          S[p,r] * fourth_mom(S, c(q,s,u,w)) * S[v,z] -
          S[p,r] * fourth_mom(S, c(q,s,u,z)) * S[v,w] -
          S[p,s] * fourth_mom(S, c(q,r,u,w)) * S[v,z] +
          S[p,s] * fourth_mom(S, c(q,r,u,z)) * S[v,w]
        ) +
        (
          S[u,w] * fourth_mom(S, c(p,r,v,z)) * S[q,s] -
          S[u,z] * fourth_mom(S, c(p,r,v,w)) * S[q,s] -
          S[u,w] * fourth_mom(S, c(p,s,v,z)) * S[q,r] +
          S[u,z] * fourth_mom(S, c(p,s,v,w)) * S[q,r]
        ) - 
        3 * (
          (S[p,r]*S[q,s] - S[p,s]*S[q,r]) * (S[u,w]*S[v,z] - S[u,z]*S[v,w])
        )
    }
  }
  return(cov)
}


test <- function(X, S, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  Y = compute_Y(X)
  cov = compute_cov_Y(X, S)
  test_stat = sqrt(n-1) * max(abs((diag(diag(cov)**(-1/2))) %*% colMeans(Y)))
  
  # sample E sets from Z~N(0,cov)
  X = mvrnorm(E, mu=rep(0,nrow(cov)), Sigma=cov)
  
  # compute critical value
  results = apply(abs((diag(diag(cov)**(-1/2))) %*% t(X)), 2, max)
  critical_values = quantile(results, probs=1-alphas)
  
  # reject or not
  is_rejected = test_stat > critical_values
  return(c(test_stat, critical_values))
}

########################################################################################################
########################################################################################################

# Sample some data
m = 10
n = 500

H = rnorm(n, mean=0, sd=1)
beta = rep(1,m)
BETA = matrix(rep(beta,n), nrow=n, byrow=T)
error_var = rep(1,m)
errors = mvrnorm(n, mu=rep(0,m), Sigma=diag(error_var))
X = BETA * H + errors
#theoretical covariance-matric 
Sigma = beta %*% t(beta) + diag(error_var)



# Different model (arbitrary covariance matrix)
A = matrix(runif(m**2)*2-1, ncol=m) 
Sigma = t(A) %*% A
X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma) # X_(i,1, ..., X_i,m are i.i.d. standard normal

# Compute sample covariance
sample_cov = matrix(0, nrow=m, ncol=m)
for (i in 1:n){
  sample_cov = sample_cov + X[i,] %*% t(X[i,])
}
sample_cov = sample_cov/n


### Sample_cov or Sigma ????
test(X, sample_cov, E=1000, alphas=0.05)




