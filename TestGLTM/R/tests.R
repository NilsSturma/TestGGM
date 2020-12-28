test_grouping <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    N = findn(nrow(X),2) # nr of estimators
    indices_U = matrix(1:N, ncol=2)
    H = calculate_H_eq(X, indices_U, ind_eq)
  } else {
    test_ineqs = TRUE
    N = findn(nrow(X),4) # nr of estimators
    indices_U = matrix(1:N, ncol=4)
    H = calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  }
  
  n = dim(H)[1]
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  
  # Mean and centering
  H_mean = Rfast::colmeans(H)
  H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  
  # Bootstrapping 
  W = bootstrap(E, H_centered)
  W[,1:p_eq] = abs(W[,1:p_eq])
  W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
  results = Rfast::rowMaxs(W_standardized, value = TRUE)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}




test_U_stat <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  
  n = dim(X)[1]  # nr of samples
  
  
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    r = 2
  } else {
    test_ineqs = TRUE
    r = 4
  }
  
  
  N = min(0.7*choose(n,r), N)
  
  # determine N_hat by Bernoulli sampling
  N_hat = rbinom(1, choose(n,r), (N / choose(n,r)))
  
  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices_U = random_combs(n,r,N_hat)
  
  # Compute matrix H
  if (test_ineqs){
    H = calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  } else {
    H = calculate_H_eq(X, indices_U, ind_eq)
  }
  
  H_mean = Rfast::colmeans(H)
  H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean)
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # equality constraints
  
  # Compute matrix G
  if (test_ineqs){
    G = calculate_G(X,L=(r-1), ind_eq, ind_ineq1, ind_ineq2)
  } else {
    G = calculate_G_eq(X, L=(r-1), ind_eq)
  }
  G_mean = Rfast::colmeans(G)
  G_centered = Rfast::transpose(Rfast::transpose(G) - G_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
  cov_G_diag = Rfast::colsums(G_centered**2) / n
  cov_diag = cov_H_diag  # (r**2) * cov_G_diag + (n/N) * cov_H_diag
  
  # Vector for standardizing
  standardizer = cov_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrap
  bootstrap_res = bootstrap_U(E, r, H_centered, G_centered, N)
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





test_U_stat_degenerate <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  
  n = dim(X)[1]  # nr of samples
  
  
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    m = 2
  } else {
    test_ineqs = TRUE
    r = 4
  }
  
  # determine N_hat by Bernoulli sampling
  N_hat = rbinom(1, choose(n,r), (N / choose(n,r)))
  
  # Choose randomly N_hat unique subsets with cardinality 4 of {1,...,n} 
  indices_U = random_combs(n,r,N_hat)
  
  # Compute matrix H
  if (test_ineqs){
    H = calculate_H(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  } else {
    H = calculate_H_eq(X, indices_U, ind_eq)
  }
  
  H_mean = Rfast::colmeans(H)
  H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean)
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # equality constraints
  
  
  # Diagonal of the sample covariance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean  ## OR N_hat??
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrap
  W = bootstrap(E, H_centered)
  W[,1:p_eq] = abs(W[,1:p_eq])
  W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
  results = Rfast::rowMaxs(W_standardized, value = TRUE)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}





test_run_over <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, B=5, E=1000, alphas=seq(0.01, 0.99, 0.01)){
  
  # Call function to calculate matrix H
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    indices_U = matrix(c(1:(nrow(X)-1),2:nrow(X)), ncol=2, byrow=FALSE)
    H = calculate_H_not_symmetric_eq(X, indices_U, ind_eq)
  } else {
    test_ineqs = TRUE
    indices_U = matrix(c(1:(nrow(X)-3),2:(nrow(X)-2),3:(nrow(X)-1),4:nrow(X)),
                       ncol=4, byrow=FALSE)
    H = calculate_H_not_symmetric(X, indices_U, ind_eq, ind_ineq1, ind_ineq2)
  }
  
  n = dim(H)[1]  # nr of samples
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # nr of equality constraints
  omega = floor(n/B)
  
  # Mean and centering
  H_mean = Rfast::colmeans(H)
  H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the batched mean estimator of the covariance of H
  cov_H_diag = rep(0, p)
  for (b in 1:omega){
    L = seq(1+(b-1)*B, b*B)
    cov_H_diag = cov_H_diag + Rfast::colsums(H_centered[L,])**2
  }
  cov_H_diag = cov_H_diag / (B*omega)
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrapping 
  W = bootstrap_m_dep(E, B, omega, H_centered)
  W[,1:p_eq] = abs(W[,1:p_eq])
  W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
  results = Rfast::rowMaxs(W_standardized, value = TRUE)
  
  # Critical values
  critical_values = quantile(results, probs=1-alphas)
  
  # Reject if test_stat > critical_value
  is_rejected = test_stat > critical_values
  return(is_rejected)
}



# test_two_step  <- function(X, ind_eq, ind_ineq1, ind_ineq2, E=1000, beta=0.001, alphas=seq(0.01, 0.99, 0.01)){
#   # Call function to calculate matrix H
#   H = calculate_H_independent(X, ind_eq, ind_ineq1, ind_ineq2)
#   
#   n = dim(H)[1]  # nr of samples
#   p = dim(H)[2]  # total nr of constraints
#   p_eq = dim(ind_eq)[1]  # nr of equality constraints
#   
#   
#   # Mean and centering
#   H_mean = Rfast::colmeans(H)
#   H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
#   
#   # Diagonal of the sample covariance of H
#   cov_H_diag = Rfast::colsums(H_centered**2) / n
#   
#   # Vector for standardizing
#   standardizer = cov_H_diag**(-1/2)
#   
#   # Test statistic
#   marginal_stats = sqrt(n) * (standardizer * H_mean)
#   test_stat =  max(abs(marginal_stats[1:p_eq]), marginal_stats[(p_eq+1):p])
#   
#   # Bootstrapping - first step (just inequalities)
#   bootstrap_res = bootstrap_independent(E, standardizer[(p_eq+1):p], H_centered[,(p_eq+1):p], 0)
#   critical_beta = quantile(bootstrap_res, probs=1-beta)
#   
#   # Selection
#   J_eq = rep(TRUE, p_eq)  # marginal_stats[1:p_eq] > -2 * critical_beta
#   J_ineq = marginal_stats[(p_eq+1):p] > -2 * critical_beta
#   p_eq = sum(J_eq)  # updated nr of equality constraints
#   p = p_eq + sum(J_ineq)  # updated total nr of constraints
#   print(p-p_eq)
#   
#   # Bootstrapping - second step
#   bootstrap_res = bootstrap_independent(E, standardizer[c(J_eq, J_ineq)], H_centered[,c(J_eq, J_ineq)], p_eq)
#   
#   # Critical values
#   critical_values = quantile(bootstrap_res, probs=(1-alphas+2*beta))
#   
#   # Reject if test_stat > critical_value
#   is_rejected = test_stat > critical_values
#   return(is_rejected)
# }




