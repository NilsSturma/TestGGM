library(MASS) #mvrnorm
library(CombMSC) #subsets




test_half_and_half <- function(X, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Testing the star tree (method: dividing X into two datasets of size n/2)
  
  n = dim(X)[1]  # Could add a test that n/2 is a natural number (or use floor(n/2))
  m = dim(X)[2]

  # Divide X into X1 and X2, both of size n/2
  X1 = X[1:(n/2),]
  X2 = X[((n/2)+1):n, ]
  
  # Compute Y_i, i=1,...,(n-1) in  a matrix
  sub_sets = subsets(m,4,1:m)
  Y = matrix(0, nrow = n/2, ncol = 2 * choose(m,4))
  for (j in 1:nrow(sub_sets)){
    p = sub_sets[j,1]
    q = sub_sets[j,2]
    r = sub_sets[j,3]
    s = sub_sets[j,4]
    Y[,j*2-1] = X1[,p] * X1[,s] * X2[,q] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
    Y[,j*2] = X1[,p] * X1[,q] * X2[,s] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
  }
  
  # Test statistic
  test_stat = max(abs( (1/sqrt((n/2))) * colSums(Y) ))  # We need absolute values here to have a two sided test
  
  # Bootstrapping 
  results = rep(0, E)
  for (i in 1:E){
    epsilons = rnorm((n/2), mean=0, sd=1)
    results[i] = max(abs( (1/sqrt((n/2))) * colSums(Y*epsilons) ))  # We need absolute values here to have a two sided tes
  }
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}




test_1_dependent <- function(X, B=3, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Testing the star tree (method similar to Dennis Leung)
  
  n = dim(X)[1]
  m = dim(X)[2]
  omega = floor((n-1)/B)
  
  # Compute Y_i, i=1,...,(n-1) and Y_mean 
  sub_sets = subsets(m,4,1:m)
  Y = matrix(0, nrow = n-1, ncol = 2 * choose(m,4))
  for (j in 1:nrow(sub_sets)){
    p = sub_sets[j,1]
    q = sub_sets[j,2]
    r = sub_sets[j,3]
    s = sub_sets[j,4]
    Y[,j*2-1] = X[1:(n-1),p] * X[1:(n-1),s] * X[2:n,q] * X[2:n,r] - X[1:(n-1),p] * X[1:(n-1),r] * X[2:n,q] * X[2:n,s]
    Y[,j*2] = X[1:(n-1),p] * X[1:(n-1),q] * X[2:n,s] * X[2:n,r] - X[1:(n-1),p] * X[1:(n-1),r] * X[2:n,q] * X[2:n,s]
  }
  
  Y_mean = colMeans(Y)
  
  # Compute the batched mean estimator for diag(cov(Y))
  cov_Y_diag = rep(0, length(Y_mean))
  Y_stand = t(t(Y) - Y_mean)  # Standardize each Y_i (Y_i - Y_mean) and save it in a matrix
  for (b in 1:omega){
    L = seq(1+(b-1)*B, b*B)
    cov_Y_diag = cov_Y_diag + colSums(Y_stand[L,])**2
  }
  cov_Y_diag = cov_Y_diag / (B*omega)
  
  # Test statistic
  test_stat = sqrt(n-1) * max(abs( diag(cov_Y_diag**(-1/2)) %*% Y_mean ))

  # Bootstrapping
  results = rep(0, E)
  for (i in 1:E){
    epsilons = rnorm(omega, mean=0, sd=1)
    sum = rep(0, length(Y_mean))
    for (b in 1:omega){
      L = seq(1+(b-1)*B, b*B)
      sum = sum + epsilons[b] * colSums(Y_stand[L,])
    }
    results[i] = max(abs( (1/sqrt(B*omega)) * diag(cov_Y_diag**(-1/2)) %*% sum ))
  }
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}




# Function to compute fourth moments (covariance matrix given)
fourth_mom = function(S, ind){
  # S: Covariance matrix
  # ind: Vector of length 4. Indices to compute the fourth moment.  
  # Example: ind=c(1,2,3,4) -> compute moment E[X_1*X_2*X_3*X_4]
  S[ind[1],ind[2]] * S[ind[3],ind[4]] + S[ind[1],ind[3]] * S[ind[2],ind[4]] + S[ind[1],ind[4]] * S[ind[2],ind[3]]
}




# Function to calculate Y: ((n-1)x(2*choose(m,4))-matrix)
# 1-dependent unbiased estimators
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




# Function to calculate covariance of Y, using the theoretical or emprical covariance matrix S
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
  for (i in 1:(2*choose(m,4))){
    p = indices[i,1]
    q = indices[i,2]
    r = indices[i,3]
    s = indices[i,4]
    for (j in i:(2*choose(m,4))){
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
      cov[j,i] = cov[i,j]
    }
  }
  return(cov)
}




test_calculate_Y <- function(X, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  m = dim(X)[2]
  n = dim(X)[1]
  
  # Sample covariance 
  S = matrix(0, nrow=m, ncol=m)
  for (i in 1:n){
    S = S + X[i,] %*% t(X[i,])
  }
  S = S/n
  
  Y = compute_Y(X)
  cov = compute_cov_Y(X, S)
  test_stat = sqrt(n-1) * max(abs((diag(diag(cov)**(-1/2))) %*% colMeans(Y)))
  
  # sample E sets from Z~N(0,cov)
  Z = mvrnorm(E, mu=rep(0,nrow(cov)), Sigma=cov)
  
  # Critical value
  results = apply(abs((diag(diag(cov)**(-1/2))) %*% t(Z)), 2, max)
  critical_values = quantile(results, probs=1-alphas)
  
  # reject or not
  is_rejected = test_stat > critical_values
  return(is_rejected)
}