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