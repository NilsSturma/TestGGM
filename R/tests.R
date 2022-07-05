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
#' @export
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
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = colSums(H_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(n) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  
  # Bootstrapping 
  W = bootstrap(E, H_centered)
  W[,1:p_eq] = abs(W[,1:p_eq])
  W_standardized = t(t(W) * standardizer)
  results = matrixStats::rowMaxs(W_standardized)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}

test_blockwise <- function(X, block_length=5, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, E=1000){

  # Nr of samples
  n = dim(X)[1]  
  
  # Check whether we test only equality constraints or not
  if (is.null(ind_ineq1) | is.null(ind_ineq2)){
    if (!is.null(ind_ineq1) | !is.null(ind_ineq2)){
      stop("ERROR - exactly one set of inequalities is missing. Cannot handle this.")
    }
    test_ineqs = FALSE
    r = 2
    p = nrow(ind_eq)
    p_eq = p
  } else {
    test_ineqs = TRUE
    r = 4
    p = nrow(ind_eq) + nrow(ind_ineq1) + nrow(ind_ineq2)
    p_eq = nrow(ind_eq)
  }

  # Get indices for subsets used in test statistic
  indices = matrix(1:findn(n,block_length), ncol=block_length)
  q = nrow(indices)
  n_block = choose(block_length,r)
  ind_large = matrix(0, nrow=(q*n_block), ncol=r)
  for (i in 1:q){
    ind_block = t(utils::combn(indices[i,], r))
    ind_large[((i-1)*n_block+1):((i-1)*n_block+n_block),] = ind_block 
  }

  # Compute estimate for each block (U-statistic for each block)
  if (test_ineqs){
    H_large = calculate_H(X, ind_large, ind_eq, ind_ineq1, ind_ineq2)
  } else {
    H_large = calculate_H_eq(X, ind_large, ind_eq)
  }
  m = matrix(0, nrow=q, ncol=nrow(ind_large))
  for (i in 1:q){
    m[i,((i-1)*n_block+1):((i-1)*n_block+n_block)] = 1
  }
  H = (1/n_block) * m %*% H_large
  
  # Mean and centering
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = colSums(H_centered**2) / q
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = sqrt(q) * H_mean
  marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrapping 
  W = bootstrap(E, H_centered)
  W[,1:p_eq] = abs(W[,1:p_eq])
  W_standardized = t(t(W) * standardizer)
  results = matrixStats::rowMaxs(W_standardized)
  
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
#' @export
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
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the batched mean estimator of the covariance of H
  cov_H_diag = rep(0, p)
  for (b in 1:omega){
    L = seq(1+(b-1)*B, b*B)
    cov_H_diag = cov_H_diag + colSums(H_centered[L,])**2
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
  W_standardized = t(t(W) * standardizer)
  results = matrixStats::rowMaxs(W_standardized)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}




#' Tests the goodness-of-fit of a Gaussian latent tree model by the maximum of a high-dimensional incomplete U-statistic
#' 
#' This function tests the goodness-of-fit of a given Gaussian latent tree model to observed data.
#' The parameter space of the model is characterized by polynomial constraints. The involved polynomials are 
#' estimable by a kernel function and the test statistic is the maximum of the corresponding incomplete U-statistic.
#'  A Gaussian multiplier bootstrap is used to estimate the limiting distribution 
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
#' @param n1 Integer, cardinality of the set \code{S_1} to compute the divide-and-conquer estimators. 
#' The subset \code{S_1} is chosen randomly from \code{\{1,...,n\}}.
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
#' @export
test_U_stat <- function(X, ind_eq, ind_ineq1=NULL, ind_ineq2=NULL, N=5000, E=1000, n1=500){
  
  
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
  
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean)
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # equality constraints
  
  # Compute matrix G on subset of samples
  X_G = X[sample(n, size=min(n1,n), replace=FALSE),]
  if (test_ineqs){
    G = calculate_G(X_G,L=(r-1), ind_eq, ind_ineq1, ind_ineq2)
  } else {
    G = calculate_G_eq(X_G, L=(r-1), ind_eq)
  }
  G_mean = colMeans(G)
  G_centered = t(t(G) - G_mean)
  
  # Diagonal of the approximate variance of H
  cov_H_diag = colSums(H_centered**2) / N_hat
  cov_G_diag = colSums(G_centered**2) / n1
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
  U_standardized = t(t(U) * standardizer)
  results = matrixStats::rowMaxs(U_standardized)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}

#########
## BSS ##
#########

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
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean) # Centering: H_i = (H_i - H_mean)

  # Diagonal of the sample covariance of H
  cov_H_diag = colSums(H_centered**2) / n

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
    W_standardized = t(t(W) * standardizer[(p_eq+1):p])
    bootstrap_res = matrixStats::rowMaxs(W_standardized)
    c_beta = as.numeric(stats::quantile(bootstrap_res, probs=1-betas[i]))


    # Calculate nuisance parameter lambda (zero for all equalities)
    lambda = H_mean[(p_eq+1):p] + (cov_H_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
    lambda = sapply(lambda, min_zero)

    # Bootstrapping - second step
    W = bootstrap(E, H_centered)
    W[,1:p_eq] = abs(W[,1:p_eq])
    W[,(p_eq+1):p] = t(t(W[,(p_eq+1):p]) + lambda * sqrt(n))
    W_standardized = t(t(W) * standardizer)
    maxima = matrixStats::rowMaxs(W_standardized)
    critical_value = as.numeric(stats::quantile(maxima, probs=1-alphas[i]+betas[i]))

    # Reject?
    results[i] = (test_stat > critical_value)
  }

  return(results)
}


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

  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean)
  p = dim(H)[2]  # total nr of constraints
  p_eq = dim(ind_eq)[1]  # equality constraints

  # Compute matrix G
  if (test_ineqs){
    G = calculate_G(X,L=(r-1), ind_eq, ind_ineq1, ind_ineq2)
  } else {
    G = calculate_G_eq(X, L=(r-1), ind_eq)
  }
  G_mean = colMeans(G)
  G_centered = t(t(G) - G_mean)

  # Diagonal of the approximate variance of H
  cov_H_diag = colSums(H_centered**2) / N_hat
  cov_G_diag = colSums(G_centered**2) / n
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
    U_standardized = t(t(U[,(p_eq+1):p]) * standardizer[(p_eq+1):p])
    bootstrap_res = matrixStats::rowMaxs(U_standardized)
    c_beta = as.numeric(stats::quantile(bootstrap_res, probs=1-betas[i]))


    # Calculate nuisance parameter lambda (zero for all equalities)
    lambda = H_mean[(p_eq+1):p] + (cov_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
    lambda = sapply(lambda, min_zero)

    # Bootstrapping - second step
    bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
    U_A = bootstrap_res[[1]]
    U_B = bootstrap_res[[2]]
    U = U_A + sqrt(n/N) * U_B
    U[,1:p_eq] = abs(U[,1:p_eq])
    U[,(p_eq+1):p] = t(t(U[,(p_eq+1):p]) + lambda * sqrt(n))
    U_standardized = t(t(U) * standardizer)
    maxima = matrixStats::rowMaxs(U_standardized)
    critical_value = as.numeric(stats::quantile(maxima, probs=1-alphas[i]+betas[i]))

    # Reject?
    results[i] = (test_stat > critical_value)
  }

  return(results)
}


#########################
## Complete Ustat Test ##
#########################

test_complete_Ustat <- function(X, ind_eq, E=1000){

    n = dim(X)[1]
    p = dim(ind_eq)[1]
    r = 2

    # Sample covariance
    S = (1/n) * t(X) %*% X

    # Complete U-statistic
    U = rep(0, p)
    for (j in 1:p){
        U[j] = S[ind_eq[j,1], ind_eq[j,2]] * S[ind_eq[j,3], ind_eq[j,4]] - 
            S[ind_eq[j,5], ind_eq[j,6]] * S[ind_eq[j,7], ind_eq[j,8]]
    }
    U = (n/(n-1)) * U  # scaling to remove bias

    # Compute matrix G
    G = calculate_G_eq(X, L=1, ind_eq)
    G_mean = colMeans(G)
    G_centered = t(t(G) - G_mean)  # note: U == G_mean

    # Diagonal of the approximate covariance matrix
    cov_diag = r**2 * (1/n) * colSums(G_centered**2)

    # Vector for standardizing
    standardizer = cov_diag**(-1/2)

    # Test statistic
    test_stat = max(abs(sqrt(n) * standardizer * U))

    # Bootstrap
    W = abs(r * bootstrap(E, G_centered))
    W_standardized = t(t(W) * standardizer)
    results = matrixStats::rowMaxs(W_standardized)
    
    # pval
    pval = (1 + sum(results >= test_stat)) / (1+E)

    return(list("PVAL"=pval, "TSTAT"=test_stat))
}
