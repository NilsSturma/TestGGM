library(CombMSC) # subsets
library(Rfast) # transpose, colmeans, rmvnorm
library(RcppHelpers)



test_half_and_half <- function(X, E=1000, alphas=seq(0.01, 0.99, 0.01), seed){
  
  # Testing the star tree (method: dividing X into two datasets of size n/2)
  
  n = dim(X)[1]  # Could add a test that n/2 is a natural number (or use floor(n/2))
  m = dim(X)[2]

  # Divide X into X1 and X2, both of size n/2
  X1 = X[1:(n/2),]
  X2 = X[((n/2)+1):n, ]
  
  # Create indices
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2*nrow(sub_sets)
  indices = matrix(0, nrow=nr_cols, 4)
  indices[1:(nr_cols/2),] = sub_sets
  indices[((nr_cols/2)+1):nr_cols,] = sub_sets[,c(1,3,2,4)]
  mode(indices) = "integer"
  
  # Compute Y
  Y = calculate_Y(indices, X1, X2) # each column is one polynom
  Y_mean = colmeans(Y)
  Y_centered = transpose(transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Diagonal of the sample covariance of Y
  cov_Y_diag = colsums(Y_centered**2) / (n/2)
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  test_stat = sqrt(n/2) * max(abs(standardizer * Y_mean))  # We need absolute values here to have a two sided test
  
  # Bootstrapping 
  results = bootstrap_independent(E, standardizer, Y_centered)
  
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
  
  # Create indices
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2*nrow(sub_sets)
  indices = matrix(0, nrow=nr_cols, 4)
  indices[1:(nr_cols/2),] = sub_sets
  indices[((nr_cols/2)+1):nr_cols,] = sub_sets[,c(1,3,2,4)]
  mode(indices) = "integer"
  
  # Compute Y
  Y = calculate_Y(indices, X1, X2) # each column is one polynom
  Y_mean = colmeans(Y)
  Y_centered = transpose(transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Compute the diagonal of the batched mean estimator cov(Y) # could do this in C++ as well (marginal effect)
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
  results = bootstrap_dependent(E, B, omega, standardizer, Y_centered)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}




test_calculate_Y <- function(X, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  m = dim(X)[2]
  n = dim(X)[1]
  X1 = X[1:(n-1),]
  X2 = X[2:n,]
  
  # Sample covariance 
  S = matrix(0, nrow=m, ncol=m)
  for (i in 1:n){
    S = S + X[i,] %*% t(X[i,])
  }
  S = S/n
  
  # create all indices
  sub_sets = subsets(m,4,1:m)
  nr_cols = 2*nrow(sub_sets)
  indices = matrix(0, nrow=nr_cols, 4)
  indices[1:(nr_cols/2),] = sub_sets
  indices[((nr_cols/2)+1):nr_cols,] = sub_sets[,c(1,3,2,4)]
  mode(indices) = "integer"
  
  # Call functions to calculate Y and its estimated covariance
  Y = calculate_Y(indices, X1, X2)
  cov = compute_cov_Y(S, indices)
  
  # Standardizer
  standardizer = diag(cov)**(-1/2)
  
  # Test statictic
  test_stat = sqrt(n-1) * max(abs(standardizer * colmeans(Y)))
  
  # Sample E sets from Z~N(0,cov)
  Z = mvrnorm(E, mu=rep(0,nrow(cov)), Sigma=cov)
  # Z is a matrix of dim=(E, nrow(cov))
  
  # Critical value
  results = apply(abs(standardizer * transpose(Z)), 2, max)
  critical_values = quantile(results, probs=1-alphas)
  
  # reject or not
  is_rejected = test_stat > critical_values
  return(is_rejected)
}