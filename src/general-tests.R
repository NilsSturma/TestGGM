test_equality_constraints = function(X, Q, not_Q, method="1-dependent", B=3, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # TODO: Add other methods
  # Remark: Instead of "Q" and "not_Q" one could give the graph itself and then call the function "findQ".
  
  n = dim(X)[1]
  omega = floor((n-1)/B)
  
  len_Q = length(Q)
  len_not_Q = length(not_Q)
  
  # Compute Y
  len_Y = length(Q) + 2*length(not_Q)
  Y = matrix(0, nrow=(n-1), ncol=len_Y)
  for (i in seq_along(Q)){
    p = Q[[i]][1]
    q = Q[[i]][2]
    r = Q[[i]][3]
    s = Q[[i]][4]
    Y[,i] = X[1:(n-1),p] * X[1:(n-1),r] * X[2:n,q] * X[2:n,s] - X[1:(n-1),p] * X[1:(n-1),s] * X[2:n,q] * X[2:n,r]
  }
  for (i in seq_along(not_Q)){
    p = not_Q[[i]][1]
    q = not_Q[[i]][2]
    r = not_Q[[i]][3]
    s = not_Q[[i]][4]
    Y[,(len_Q+i)] = X[1:(n-1),p] * X[1:(n-1),r] * X[2:n,q] * X[2:n,s] - X[1:(n-1),p] * X[1:(n-1),s] * X[2:n,q] * X[2:n,r]
    Y[,(len_Q + len_not_Q +i)] = X[1:(n-1),p] * X[1:(n-1),q] * X[2:n,r] * X[2:n,s] - X[1:(n-1),p] * X[1:(n-1),s] * X[2:n,q] * X[2:n,r]
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
  return(c(test_stat, critical_values))
}