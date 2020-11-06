library(CombMSC) #subsets
library(Rfast) #Rnorm, transpose, colmeans, rmvnorm



test_half_and_half <- function(X, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Testing the star tree (method: dividing X into two datasets of size n/2)
  
  n = dim(X)[1]  # Could add a test that n/2 is a natural number (or use floor(n/2))
  m = dim(X)[2]

  # Divide X into X1 and X2, both of size n/2
  X1 = X[1:(n/2),]
  X2 = X[((n/2)+1):n, ]
  
  # Compute Y_i, i=1,...,(n/2) in  a matrix
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2 * nrow(sub_sets)
  Y = matrix(0, nrow = n/2, ncol = nr_cols)
  for (j in 1:(nr_cols/2)){
    p = sub_sets[j,1]
    q = sub_sets[j,2]
    r = sub_sets[j,3]
    s = sub_sets[j,4]
    Y[,j*2-1] = X1[,p] * X1[,s] * X2[,q] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
    Y[,j*2] = X1[,p] * X1[,q] * X2[,s] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
  }
  
  Y_mean = colmeans(Y)
  Y_centered = transpose(transpose(Y) - Y_mean) # Center each Y_i = (Y_i - Y_mean) and save it in a matrix
  
  # Diagonal of the sample covariance of Y
  cov_Y_diag = colsums(Y_centered**2) / (n/2)
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  test_stat = sqrt(n/2) * max(abs(standardizer * Y_mean))  # We need absolute values here to have a two sided test
  
  # Bootstrapping 
  results = rep(0, E)
  for (i in 1:E){
    epsilons = Rnorm((n/2), m=0, s=1)
    results[i] = (1/sqrt((n/2))) * max(abs(standardizer * colsums(Y_centered*epsilons)))  # We need absolute values here to have a two sided tes
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
  X1 = X[1:(n-1),]
  X2 = X[2:n,]
  
  # Compute Y_i, i=1,...,(n-1) and Y_mean 
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2 * nrow(sub_sets)
  Y = matrix(0, nrow = n-1, ncol = nr_cols)
  for (j in 1:(nr_cols/2)){
    p = sub_sets[j,1]
    q = sub_sets[j,2]
    r = sub_sets[j,3]
    s = sub_sets[j,4]
    Y[,j*2-1] = X1[,p] * X1[,s] * X2[,q] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
    Y[,j*2] =  X1[,p] * X1[,q] * X2[,s] * X2[,r] - X1[,p] * X1[,r] * X2[,q] * X2[,s]
  }
  
  Y_mean = colmeans(Y)
  Y_centered = transpose(transpose(Y) - Y_mean)  # Center each Y_i = (Y_i - Y_mean) and save it in a matrix
  
  # Compute the diagonal of the batched mean estimator cov(Y)
  cov_Y_diag = rep(0, length(Y_mean))
  for (b in 1:omega){
    L = seq(1+(b-1)*B, b*B)
    cov_Y_diag = cov_Y_diag + colsums(Y_centered[L,])**2
  }
  cov_Y_diag = cov_Y_diag / (B*omega)
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  test_stat = sqrt(n-1) * max(abs(standardizer * Y_mean))

  # Bootstrapping
  results = rep(0, E)
  for (i in 1:E){
    epsilons = Rnorm(omega, m=0, s=1)
    sum = rep(0, length(Y_mean))
    for (b in 1:omega){
      L = seq(1+(b-1)*B, b*B)
      sum = sum + epsilons[b] * colsums(Y_centered[L,])
    }
    results[i] = (1/sqrt(B*omega)) * max(abs(standardizer * sum))
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
compute_Y = function(X, indices){
  
  m = dim(X)[2]
  n = dim(X)[1]
  nr_cols = nrow(indices)
  
  # Compute Y (unbiased estimate of the tetrads)
  Y = matrix(0, nrow = n-1, ncol=nr_cols)
  
  for (i in 1:nr_cols){
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
compute_cov_Y = function(S, indices){
  
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
  
  # create all indices in right order !!!different order than previously!!!
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2*nrow(sub_sets)
  indices = matrix(0, nrow=nr_cols, 4)
  indices[1:(nr_cols/2),] = sub_sets
  indices[((nr_cols/2)+1):(2*(nr_cols/2)),] = sub_sets[,c(1,3,2,4)]
  
  # call functions
  Y = compute_Y(X, indices)
  cov = compute_cov_Y(S, indices)
  
  # Standardizer
  standardizer = diag(cov)**(-1/2)
  
  # Test statictic
  test_stat = sqrt(n-1) * max(abs(standardizer * colmeans(Y)))
  
  # sample E sets from Z~N(0,cov)
  Z = rmvnorm(E, mu=rep(0,nrow(cov)), sigma=cov)
  # Z is a matrix of dim=(E, nrow(cov))
  
  # Critical value
  results = apply(abs(standardizer * transpose(Z)), 2, max)
  critical_values = quantile(results, probs=1-alphas)
  
  # reject or not
  is_rejected = test_stat > critical_values
  return(is_rejected)
}