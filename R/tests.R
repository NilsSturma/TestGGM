#' Tests the goodness-of-fit of a Gaussian latent tree model by the maximum of a high-dimensional independent sum
#' 
#' This function tests the goodness-of-fit of a given Gaussian latent tree model to observed data.
#' The parameter space is characterized by polynomial constraints. The involved polynomials are 
#' estimated by grouping the data into independent subsets. Each group is used to form an 
#' unbiased estimate of all polynomials. The test statistic 
#' is the maximum of the average of the independent studentized estimates.
#' A Gaussian multiplier bootstrap is used to estimate the limiting distribution of the 
#' test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. The number of columns has to be equal to the number of 
#' leaves of the postulated tree (i.e. number of observed variables). Each row corresponds to one sample.
#' @param ind_eq Representation of the equality constraints that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}.
#' @param ind_ineq1 Representation of the inequality constraints in three variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param ind_ineq2 Representation of the inequality constraints in four variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param E Integer, number of bootstrap iterations.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Create tree
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' 
#' # Sample data from tree
#' igraph::V(tree)$var = rep(1,8)
#' igraph::E(tree)$corr = rep(0.7,7)
#' X = sample_from_tree(tree, m=5, n=500)
#' 
#' # Determine the representation of the polynomials that have to be tested
#' res = collect_indices(tree, m=5)
#' ind_eq = res$ind_eq
#' ind_ineq1 = res$ind_ineq1
#' ind_ineq2 = res$ind_ineq2
#' 
#' # Apply the test
#' test_grouping(X, ind_eq, ind_ineq1, ind_ineq2)
#' 
#' @references
#' TO BE WRITTEN
test_grouping <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, E=1000){
  
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    N = findn(nrow(X),2)
    indices = matrix(1:N, ncol=2)
    H = calculate_H_eq(X, indices, ind_eq)
  } else {
    test_ineqs = TRUE
    N = findn(nrow(X),4)
    indices = matrix(1:N, ncol=4)
    H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
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
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}





#' Tests the goodness-of-fit of a Gaussian latent tree model by the maximum of a high-dimensional m-dependent sum
#' 
#' This function tests the goodness-of-fit of a given Gaussian latent tree model to observed data.
#' The parameter space is characterized by polynomial constraints. The involved polynomials are 
#' estimated by considering overlapping subsets of the data 
#' (\code{\{X_1, X_2, X_3, X_4\}, \{X_2, X_3, X_4, X_5\}, ...}). Each subset is used to form an 
#' unbiased estimate of all polynomials. The test statistic 
#' is the maximum of the average of the m-dependent studentized estimates.
#' A Gaussian multiplier bootstrap is used to estimate the limiting distribution of the 
#' test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. The number of columns has to be equal to the number of 
#' leaves of the postulated tree (i.e. number of observed variables). Each row corresponds to one sample.
#' @param ind_eq Representation of the equality constraints that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}.
#' @param ind_ineq1 Representation of the inequality constraints in three variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param ind_ineq2 Representation of the inequality constraints in four variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param B Integer, batch size for the estimate of the covariance matrix.
#' @param E Integer, number of bootstrap iterations.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Create tree
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' 
#' # Sample data from tree
#' igraph::V(tree)$var = rep(1,8)
#' igraph::E(tree)$corr = rep(0.7,7)
#' X = sample_from_tree(tree, m=5, n=500)
#' 
#' # Determine the representation of the polynomials that have to be tested
#' res = collect_indices(tree, m=5)
#' ind_eq = res$ind_eq
#' ind_ineq1 = res$ind_ineq1
#' ind_ineq2 = res$ind_ineq2
#' 
#' # Apply the test
#' test_run_over(X, ind_eq, ind_ineq1, ind_ineq2)
#' 
#' @references
#' TO BE WRITTEN
test_run_over <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, B=5, E=1000){
  
  # Call function to calculate matrix H
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    indices = matrix(c(1:(nrow(X)-1),2:nrow(X)), ncol=2, byrow=FALSE)
    H = calculate_H_not_symmetric_eq(X, indices, ind_eq)
  } else {
    test_ineqs = TRUE
    indices = matrix(c(1:(nrow(X)-3),2:(nrow(X)-2),3:(nrow(X)-1),4:nrow(X)),
                       ncol=4, byrow=FALSE)
    H = calculate_H_not_symmetric(X, indices, ind_eq, ind_ineq1, ind_ineq2)
  }
  
  n = dim(H)[1]  # nr of estimates
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
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}




#' Tests the goodness-of-fit of a Gaussian latent tree model by the maximum of a high-dimensional U-statistic
#' 
#' This function tests the goodness-of-fit of a given Gaussian latent tree model to observed data.
#' The parameter space of the model is characterized by polynomial constraints. The involved polynomials are 
#' estimated by considering subsets of the data. The number of subsets as well as the subsets itself 
#' are chosen randomly. Each subset is used to form an unbiased estimate of the polynomials. 
#' The test statistic is the maximum of the U-statistic formed by 
#' the studentized estimates. A Gaussian multiplier bootstrap is used to estimate the limiting distribution 
#' of the  test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. The number of columns has to be equal to the number of 
#' leaves of the postulated tree (i.e. number of observed variables). Each row corresponds to one sample.
#' @param ind_eq Representation of the equality constraints that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}.
#' @param ind_ineq1 Representation of the inequality constraints in three variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param ind_ineq2 Representation of the inequality constraints in four variables that have to be tested. 
#' Create this object using the function \code{\link{collect_indices}}. If \code{NULL} inequality constraints are not tested.
#' @param N Integer, computational budget parameter.
#' @param E Integer, number of bootstrap iterations.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Create tree
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' 
#' # Sample data from tree
#' igraph::V(tree)$var = rep(1,8)
#' igraph::E(tree)$corr = rep(0.7,7)
#' X = sample_from_tree(tree, m=5, n=500)
#' 
#' # Determine the representation of the polynomials that have to be tested
#' res = collect_indices(tree, m=5)
#' ind_eq = res$ind_eq
#' ind_ineq1 = res$ind_ineq1
#' ind_ineq2 = res$ind_ineq2
#' 
#' # Apply the test
#' test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2)
#' 
#' @references
#' TO BE WRITTEN
test_U_stat <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000){
  
  
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
  N_hat = stats::rbinom(1, choose(n,r), (N / choose(n,r)))
  
  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices = matrix(unlist(random_combs(n,r,N_hat)[[1]]), ncol = r, byrow = TRUE)
  
  # Compute matrix H
  if (test_ineqs){
    H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
  } else {
    H = calculate_H_eq(X, indices, ind_eq)
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
  
  # Diagonal of the approximate variance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
  cov_G_diag = Rfast::colsums(G_centered**2) / n
  #cov_diag = cov_H_diag
  cov_diag = r**2 * cov_G_diag + (n/N) * cov_H_diag
  
  # Vector for standardizing
  standardizer = cov_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrap
  bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
  U_A = bootstrap_res[[1]]
  U_B = bootstrap_res[[2]]
  U = U_A + sqrt(n/N) * U_B
  U[,1:p_eq] = abs(U[,1:p_eq])
  U_standardized = Rfast::transpose(Rfast::transpose(U) * standardizer)
  results = Rfast::rowMaxs(U_standardized, value = TRUE)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}


# helper function
min_zero <- function(x){min(x,0)}


test_grouping_BSS  <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, E=1000, alphas=0.05, betas=0.005){

  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    N = findn(nrow(X),2)
    indices = matrix(1:N, ncol=2)
    H = calculate_H_eq(X, indices, ind_eq)
  } else {
    test_ineqs = TRUE
    N = findn(nrow(X),4)
    indices = matrix(1:N, ncol=4)
    H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
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
  
  results = rep(FALSE, length(alphas))
  for (i in 1:length(alphas)){
    
    # Bootstrapping - first step  (only inequalities)
    W = bootstrap(E, H_centered[,(p_eq+1):p])
  
    # Calculate c_beta
    W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer[(p_eq+1):p])
    bootstrap_res = Rfast::rowMaxs(W_standardized, value = TRUE)
    c_beta = as.numeric(quantile(bootstrap_res, probs=1-betas[i]))


    # Calculate nuisance parameter lambda (zero for all equalities)
    lambda = H_mean[(p_eq+1):p] + (cov_H_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
    lambda = sapply(lambda, min_zero)

    # Bootstrapping - second step
    W = bootstrap(E, H_centered)
    W[,1:p_eq] = abs(W[,1:p_eq])
    W[,(p_eq+1):p] = Rfast::transpose(Rfast::transpose(W[,(p_eq+1):p]) + lambda * sqrt(n))
    W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
    maxima = Rfast::rowMaxs(W_standardized, value = TRUE)
    critical_value = as.numeric(quantile(maxima, probs=1-alphas[i]+betas[i]))

    # Reject?
    results[i] = (test_stat > critical_value)
  }

  return(results)
}


# test_grouping_BSS  <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, E=1000, alpha=0.05, beta=0.005){
# 
#   if (is.null(ind_ineq1) | is.null(ind_ineq2)){
#     if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
#       stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
#     }
#     test_ineqs = FALSE
#     N = findn(nrow(X),2)
#     indices = matrix(1:N, ncol=2)
#     H = calculate_H_eq(X, indices, ind_eq)
#   } else {
#     test_ineqs = TRUE
#     N = findn(nrow(X),4)
#     indices = matrix(1:N, ncol=4)
#     H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
#   }
# 
#   n = dim(H)[1]
#   p = dim(H)[2]  # total nr of constraints
#   p_eq = dim(ind_eq)[1]  # nr of equality constraints
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
#   marginal_stats = sqrt(n) * H_mean
#   marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
#   test_stat =  max(standardizer * marginal_stats)
# 
#   # Bootstrapping - first step (only inequalities)
#   W = bootstrap(E, H_centered[,(p_eq+1):p])
#   W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer[(p_eq+1):p])
#   bootstrap_res = Rfast::rowMaxs(W_standardized, value = TRUE)
#   c_beta= as.numeric(quantile(bootstrap_res, probs=1-beta))
# 
#   # Calculate nuisance parameter lambda (zero for all equalities)
#   lambda = H_mean[(p_eq+1):p] + (cov_H_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
#   lambda = sapply(lambda, min_zero)
# 
#   # Bootstrapping - second step
#   W = bootstrap(E, H_centered)
#   W[,1:p_eq] = abs(W[,1:p_eq])
#   W[,(p_eq+1):p] = Rfast::transpose(Rfast::transpose(W[,(p_eq+1):p]) + lambda * sqrt(n))
#   W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
#   maxima = Rfast::rowMaxs(W_standardized, value = TRUE)
#   critical_value = as.numeric(quantile(maxima, probs=1-alpha+beta))
# 
#   # Reject?
#   result = (test_stat > critical_value)
# 
#   return(result)
# }





test_U_stat_BSS <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000, alphas=0.05, betas=0.005){


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
  N_hat = stats::rbinom(1, choose(n,r), (N / choose(n,r)))

  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices = matrix(unlist(random_combs(n,r,N_hat)[[1]]), ncol = r, byrow = TRUE)

  # Compute matrix H
  if (test_ineqs){
    H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
  } else {
    H = calculate_H_eq(X, indices, ind_eq)
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

  # Diagonal of the approximate variance of H
  cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
  cov_G_diag = Rfast::colsums(G_centered**2) / n
  cov_diag = r**2 * cov_G_diag + (n/N) * cov_H_diag

  # Vector for standardizing
  standardizer = cov_diag**(-1/2)

  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  results = rep(FALSE, length(alphas))
  for (i in 1:length(alphas)){
    
    # Bootstrap - first step
    bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
    U_A = bootstrap_res[[1]]
    U_B = bootstrap_res[[2]]
    U = U_A + sqrt(n/N) * U_B
  
    # Calculate c_beta
    U_standardized = Rfast::transpose(Rfast::transpose(U[,(p_eq+1):p]) * standardizer[(p_eq+1):p])
    bootstrap_res = Rfast::rowMaxs(U_standardized, value = TRUE)
    c_beta = as.numeric(quantile(bootstrap_res, probs=1-betas[i]))


    # Calculate nuisance parameter lambda (zero for all equalities)
    lambda = H_mean[(p_eq+1):p] + (cov_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
    lambda = sapply(lambda, min_zero)

    # Bootstrapping - second step
    bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
    U_A = bootstrap_res[[1]]
    U_B = bootstrap_res[[2]]
    U = U_A + sqrt(n/N) * U_B
    U[,1:p_eq] = abs(U[,1:p_eq])
    U[,(p_eq+1):p] = Rfast::transpose(Rfast::transpose(U[,(p_eq+1):p]) + lambda * sqrt(n))
    U_standardized = Rfast::transpose(Rfast::transpose(U) * standardizer)
    maxima = Rfast::rowMaxs(U_standardized, value = TRUE)
    critical_value = as.numeric(quantile(maxima, probs=1-alphas[i]+betas[i]))

    # Reject?
    results[i] = (test_stat > critical_value)
  }

  return(results)
}


# test_U_stat_BSS <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000, alpha=0.05, beta=0.005){
# 
# 
#   n = dim(X)[1]  # nr of samples
# 
# 
#   if (is.null(ind_ineq1) | is.null(ind_ineq2)){
#     if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
#       stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
#     }
#     test_ineqs = FALSE
#     r = 2
#   } else {
#     test_ineqs = TRUE
#     r = 4
#   }
# 
# 
#   N = min(0.7*choose(n,r), N)
# 
#   # determine N_hat by Bernoulli sampling
#   N_hat = stats::rbinom(1, choose(n,r), (N / choose(n,r)))
# 
#   # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
#   indices = matrix(unlist(random_combs(n,r,N_hat)[[1]]), ncol = r, byrow = TRUE)
# 
#   # Compute matrix H
#   if (test_ineqs){
#     H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)
#   } else {
#     H = calculate_H_eq(X, indices, ind_eq)
#   }
# 
#   H_mean = Rfast::colmeans(H)
#   H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean)
#   p = dim(H)[2]  # total nr of constraints
#   p_eq = dim(ind_eq)[1]  # equality constraints
# 
#   # Compute matrix G
#   if (test_ineqs){
#     G = calculate_G(X,L=(r-1), ind_eq, ind_ineq1, ind_ineq2)
#   } else {
#     G = calculate_G_eq(X, L=(r-1), ind_eq)
#   }
#   G_mean = Rfast::colmeans(G)
#   G_centered = Rfast::transpose(Rfast::transpose(G) - G_mean)
# 
#   # Diagonal of the approximate variance of H
#   cov_H_diag = Rfast::colsums(H_centered**2) / N_hat
#   cov_G_diag = Rfast::colsums(G_centered**2) / n
#   cov_diag = r**2 * cov_G_diag + (n/N) * cov_H_diag
# 
#   # Vector for standardizing
#   standardizer = cov_diag**(-1/2)
# 
#   # Test statistic
#   marginal_stats = sqrt(n) * H_mean
#   marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
#   test_stat =  max(standardizer * marginal_stats)
# 
#   # Bootstrap # first step
#   bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
#   U_A = bootstrap_res[[1]]
#   U_B = bootstrap_res[[2]]
#   U = U_A + sqrt(n/N) * U_B
# 
# 
#   # Calculate c_beta
#   U_standardized = Rfast::transpose(Rfast::transpose(U[,(p_eq+1):p]) * standardizer[(p_eq+1):p])
#   bootstrap_res = Rfast::rowMaxs(U_standardized, value = TRUE)
#   c_beta = as.numeric(quantile(bootstrap_res, probs=1-beta))
# 
# 
#   # Calculate nuisance parameter lambda (zero for all equalities)
#   lambda = H_mean[(p_eq+1):p] + (cov_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
#   lambda = sapply(lambda, min_zero)
# 
#   # Bootstrapping - second step
#   bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
#   U_A = bootstrap_res[[1]]
#   U_B = bootstrap_res[[2]]
#   U = U_A + sqrt(n/N) * U_B
#   U[,1:p_eq] = abs(U[,1:p_eq])
#   U[,(p_eq+1):p] = Rfast::transpose(Rfast::transpose(U[,(p_eq+1):p]) + lambda * sqrt(n))
#   U_standardized = Rfast::transpose(Rfast::transpose(U) * standardizer)
#   maxima = Rfast::rowMaxs(U_standardized, value = TRUE)
#   critical_value = as.numeric(quantile(maxima, probs=1-alpha+beta))
# 
#   # Reject?
#   result = (test_stat > critical_value)
# 
#   return(result)
# }









# ############# grouping_compute_cov ##############
# # testing only equalities
# # estimators h NOT symmetrized
# 
# test_grouping_compute_cov <- function(X, ind_eq, E=1000){
#   
#   # Call function to calculate matrix H
#   N = findn(nrow(X),2)
#   indices = matrix(1:N, ncol=2)
#   H = calculate_H_not_symmetric_eq(X, indices, ind_eq)
#   
#   # Mean and centering
#   H_mean = Rfast::colmeans(H)
#   H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
#   
#   # Nr of samples
#   n = dim(X)[1]
#   
#   # Sample covariance
#   S = (1/n)*t(X)%*%X
#   
#   # Estimate limiting covariance matrix
#   cov = cov_grouping(S, ind_eq)
#   
#   # Vector for standardizing
#   standardizer = diag(cov)**(-1/2)
#   
#   # Test statictic
#   test_stat = sqrt(n/2) * max(abs(standardizer * H_mean))
#   
#   # Sample E sets from Z~N(0,cov)
#   Z = Rfast::rmvnorm(E, mu=rep(0,nrow(cov)), sigma=cov)
#   
#   # Critical value
#   results = apply(abs(standardizer * Rfast::transpose(Z)), 2, max)
#   
#   # pval
#   pval = (1 + sum(results >= test_stat)) / (1+E)
#   
#   return(list("PVAL"=pval, "TSTAT"=test_stat))
# }
# 
# 
# 
# 
# 
# 
# 
# ############# run_over_compute_cov ##############
# # testing only equalities
# # estimators h NOT symmetrized
# 
# test_run_over_compute_cov <- function(X, ind_eq, E=1000){
#   
#   # Call function to calculate matrix H
#   indices = matrix(c(1:(nrow(X)-1),2:nrow(X)), ncol=2, byrow=FALSE)
#   H = calculate_H_not_symmetric_eq(X, indices, ind_eq)
#   
#   # Mean and centering
#   H_mean = Rfast::colmeans(H)
#   H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)
#   
#   # Nr of samples
#   n = dim(X)[1]
#   
#   # Sample covariance
#   S = (1/n)*t(X)%*%X
#   
#   # Estimate limiting covariance matrix
#   cov = cov_run_over(S, ind_eq)
#   
#   # Vector for standardizing
#   standardizer = diag(cov)**(-1/2)
#   
#   # Test statictic
#   test_stat = sqrt(n-1) * max(abs(standardizer * H_mean))
#   
#   # Sample E sets from Z~N(0,cov)
#   Z = Rfast::rmvnorm(E, mu=rep(0,nrow(cov)), sigma=cov)
#   
#   # Critical value
#   results = apply(abs(standardizer * Rfast::transpose(Z)), 2, max)
#   
#   # pval
#   pval = (1 + sum(results >= test_stat)) / (1+E)
#   
#   return(list("PVAL"=pval, "TSTAT"=test_stat))
# }





