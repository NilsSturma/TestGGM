test_independent <- function(X, ind_eq, ind_ineq1, ind_ineq2, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Call function to calculate matrix Y
  Y = calculate_Y_independent(X, ind_eq, ind_ineq1, ind_ineq2)
  
  n = dim(Y)[1]  # nr of samples
  p = dim(Y)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  
  # Mean and centering
  Y_mean = Rfast::colmeans(Y)
  Y_centered = Rfast::transpose(Rfast::transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Diagonal of the sample covariance of Y
  cov_Y_diag = Rfast::colsums(Y_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * standardizer * Y_mean
  test_stat =  max(abs(marginal_stats[1:p_eq]), marginal_stats[(p_eq+1):p])
  
  # Bootstrapping 
  bootstrap_res = bootstrap_independent(E, standardizer, Y_centered, p_eq)
  
  # Critical values
  critical_values = quantile(bootstrap_res, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}






test_two_step  <- function(X, ind_eq, ind_ineq1, ind_ineq2, E=1000, beta=0.001, alphas=seq(0.01, 0.99, 0.01)){
  # Call function to calculate matrix Y
  Y = calculate_Y_independent(X, ind_eq, ind_ineq1, ind_ineq2)
  
  n = dim(Y)[1]  # nr of samples
  p = dim(Y)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  
  
  # Mean and centering
  Y_mean = Rfast::colmeans(Y)
  Y_centered = Rfast::transpose(Rfast::transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Diagonal of the sample covariance of Y
  cov_Y_diag = Rfast::colsums(Y_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * (standardizer * Y_mean)
  test_stat =  max(abs(marginal_stats[1:p_eq]), marginal_stats[(p_eq+1):p])
  
  # Bootstrapping - first step (just inequalities)
  bootstrap_res = bootstrap_independent(E, standardizer[(p_eq+1):p], Y_centered[,(p_eq+1):p], 0)
  critical_beta = quantile(bootstrap_res, probs=1-beta)
  
  # Selection
  J_eq = rep(TRUE, p_eq)  # marginal_stats[1:p_eq] > -2 * critical_beta
  J_ineq = marginal_stats[(p_eq+1):p] > -2 * critical_beta
  p_eq = sum(J_eq)  # updated nr of equality constraints
  p = p_eq + sum(J_ineq)  # updated total nr of constraints
  print(p-p_eq)
  
  # Bootstrapping - second step
  bootstrap_res = bootstrap_independent(E, standardizer[c(J_eq, J_ineq)], Y_centered[,c(J_eq, J_ineq)], p_eq)
  
  # Critical values
  critical_values = quantile(bootstrap_res, probs=(1-alphas+2*beta))
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}


test_symmetric <- function(X, ind_eq, ind_ineq1, ind_ineq2, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # create permutations
  perm = permutations(4)
  
  # Call function to calculate matrix Y
  Y = calculate_Y_symmetric(X, ind_eq, ind_ineq1, ind_ineq2, perm)
  # Alternative (equal)
  # indices_U = matrix(1:X.nrow(), ncol=4) # could do byrow=TRUE
  # Y = calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  
  n = dim(Y)[1]  # nr of samples
  p = dim(Y)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  
  # Mean and centering
  Y_mean = Rfast::colmeans(Y)
  Y_centered = Rfast::transpose(Rfast::transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Diagonal of the sample covariance of Y
  cov_Y_diag = Rfast::colsums(Y_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * standardizer * Y_mean
  test_stat =  max(abs(marginal_stats[1:p_eq]), marginal_stats[(p_eq+1):p])
  
  # Bootstrapping 
  bootstrap_res = bootstrap_independent(E, standardizer, Y_centered, p_eq)
  
  # Critical values
  critical_values = quantile(bootstrap_res, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}




test_U_stat <- function(X, ind_eq, ind_ineq1, ind_ineq2, N=NULL, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  
  n = dim(X)[1]  # nr of samples
  if (is.null(N)){
    N = 2*n
  }
  
  # determine N_hat by Bernoulli sampling
  N_hat = rbinom(1, choose(n,4), (N / choose(n,4)))
  
  # Choose randomly N_hat subsets with cardinality 4 of {1,...,n}
  indices_U = matrix(ncol=4, nrow=N_hat) 
  for (i in 1:N_hat){
    indices_U[i,] = sort(sample(1:n, 4, replace=FALSE), decreasing=FALSE)
  } # very unlikely that we get the same indices twice
  
  # Compute matrix H
  H = calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # equality constraints
  
  # Mean and centering
  H_mean = Rfast::colmeans(H)
  H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Compute matrix G
  G = calculate_G(X,L=3, ind_eq, ind_ineq1, ind_ineq2)
  G_mean = Rfast::colmeans(G)
  G_centered = Rfast::transpose(Rfast::transpose(G) - G_mean)
  
  # Bootstrap
  bootstrap_res = bootstrap_U(E, H_centered, G_centered, N)
  U_A = bootstrap_res[[1]]
  U_B = bootstrap_res[[2]]
  U = U_A + sqrt(n/N) * U_B
  U[,1:p_eq] = abs(U[,1:p_eq])
  U_standardized = Rfast::transpose(Rfast::transpose(U) * standardizer)
  results = Rfast::rowMaxs(U_standardized, value = TRUE)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}









test_m_dep <- function(X, ind_eq, ind_ineq1, ind_ineq2, B=5, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Call function to calculate matrix Y
  Y = calculate_Y_m_dep(X, ind_eq, ind_ineq1, ind_ineq2)
  
  n = dim(Y)[1]  # nr of samples
  p = dim(Y)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  omega = floor(n/B)
  
  # Mean and centering
  Y_mean = Rfast::colmeans(Y)
  Y_centered = Rfast::transpose(Rfast::transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)
  
  # Diagonal of the batched mean estimator of the covariance of Y
  cov_Y_diag = rep(0, p)
  for (b in 1:omega){
    L = seq(1+(b-1)*B, b*B)
    cov_Y_diag = cov_Y_diag + Rfast::colsums(Y_centered[L,])**2
  }
  cov_Y_diag = cov_Y_diag / (B*omega)
  
  
  # Vector for standardizing
  standardizer = cov_Y_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * (standardizer * Y_mean)
  test_stat =  max(abs(marginal_stats[1:p_eq]), marginal_stats[(p_eq+1):p])
  
  # Bootstrapping 
  results = bootstrap_m_dep(E, B, omega, standardizer, Y_centered, p_eq)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}

