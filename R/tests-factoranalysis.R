#' Randomly choosing a set of minors
#' 
#' Randomly samples indices of minors without replacement.
#' 
#' @param m Integer, number of observed nodes (and dimension of (m x m)-matrix).
#' @param factors Integer, number of hidden factors.
#' @param nr Integer, number of (\code{factors+1 x factors+1}) minors to sample.
#' 
#' @return Matrix with \code{nr} rows and \code{2*factors+2} columns. Each row corresponds to one minor.
#' The first \code{factors+1} columns correspond to the row-indices of the minor, 
#' the last \code{factors+1} columns correspond to the column-indices of the minor.
#' 
#' @examples
#' random_minors(20,2,10000)
random_minors <- function(m, factors, nr){
  
  if ((2*(factors+1)) > m){
    stop("There is no (factors+1 x factors+1)-minor in a (m x m)-matrix.")
  }
  
  if (choose(m,(2*(factors+1))) < nr){
    stop("nr is larger than the number of possible (factors+1 x factors+1)-minors of a (m x m)-matrix.")
  }
  
  # choose sets of indices
  A = matrix(unlist(random_combs(m,(2*(factors+1)),nr)[[1]]), 
             ncol = (2*(factors+1)), byrow = TRUE)
  
  # shuffle each row
  A = t(apply(A, 1, sample))
  
  A = cbind(t(apply(A[,1:(factors+1)],1,sort)), t(apply(A[,(factors+2):(2*(factors+1))],1,sort)))
  
  return(A)
}


#' Tests the (3x3)-minor constraints of a 2-factor model by the maximum of a high-dimensional independent sum
#' 
#' This function tests the (3x3)-minor constraints of a given 2-factor model to observed data.
#' The minors are estimated by grouping the data into independent subsets. 
#' Each group is used to form an unbiased estimate of all minors. 
#' The test statistic is the maximum of the average of the independent studentized estimates. 
#' A Gaussian multiplier bootstrap procedure is used to estimate the limiting distribution of 
#' the test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. The number of columns corresponds to the number of observed variables. 
#' Each row corresponds to one sample.
#' @param nr_minors Integer, number of randomly chosen minors that are tested. 
#' @param E Integer, number of bootstrap iterations.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Covariance matrix corresponding to the two-factor analysis model
#' m=20
#' Gamma = matrix(stats::rnorm(2*m),m,2)
#' Psi = diag(rep(1,m))
#' cov = Psi + Gamma %*% t(Gamma)
#' 
#' # Sample data from the two-factor analysis model
#' X = MASS::mvrnorm(500, mu=rep(0,nrow(cov)), Sigma=cov)
#' 
#' # Apply the test
#' test_indep_factors(X, 10000)
#' 
#' @references 
#' TO BE WRITTEN
test_indep_factors <- function(X, nr_minors, E=1000){
  
  m = dim(X)[2] # nr of observed variables
  
  # split dataset
  N = findn(nrow(X),3)
  indices = matrix(1:N, ncol=3)
  
  # randomly choose minors to test
  ind_minors = random_minors(m,2,nr_minors)
  
  # Calculate matrix H of estimates
  H = H_factors(X, indices, ind_minors)
  n = dim(H)[1] # number of estimates
  p = dim(H)[2]  # nr of constraints
  
  
  # Mean and centering
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean) # Centering: H_i = (H_i - H_mean)
  
  # Diagonal of the sample covariance of H
  cov_H_diag = colSums(H_centered**2) / n
  
  # Vector for standardizing
  standardizer = cov_H_diag**(-1/2)
  
  # Test statistic
  marginal_stats = abs(sqrt(n) * H_mean)
  test_stat =  max(standardizer * marginal_stats)
  
  
  # Bootstrapping 
  W = abs(bootstrap(E, H_centered))
  W_standardized = t(t(W) * standardizer)
  results = matrixStats::rowMaxs(W_standardized)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}


#' Tests the (3x3)-minor constraints of a 2-factor model by the maximum of a high-dimensional U-statistic
#' 
#' This function tests the(3x3)-minor constraints of a given 2-factor model to observed data.
#' The minors are estimated by considering subsets of the data. The number of subsets as well as
#' the subsets itself are chosen randomly. Each subset is used to form an unbiased estimate of all minors. 
#' The test statistic is the maximum of the U-statistic formed by the studentized estimates. 
#' A Gaussian multiplier bootstrap procedure is used to estimate the limiting distribution of 
#' the test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. The number of columns corresponds to the number of observed variables. 
#' Each row corresponds to one sample.
#' @param nr_minors Integer, number of randomly chosen minors that are tested. 
#' @param N Integer, computational budget parameter.
#' @param E Integer, number of bootstrap iterations.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Covariance matrix corresponding to the two-factor analysis model
#' m=20
#' Gamma = matrix(stats::rnorm(2*m),m,2)
#' Psi = diag(rep(1,m))
#' cov = Psi + Gamma %*% t(Gamma)
#' 
#' # Sample data from the two-factor analysis model
#' X = MASS::mvrnorm(500, mu=rep(0,nrow(cov)), Sigma=cov)
#' 
#' # Apply the test
#' test_U_stat_factors(X, 10000)
#' 
#' @references 
#' TO BE WRITTEN
test_U_stat_factors <- function(X, nr_minors, N=5000, E=1000){
  
  
  n = dim(X)[1]  # nr of samples
  m = dim(X)[2] # nr of observed variables
  r = 3 # order of U-statistic
  
  N = min(0.7*choose(n,r), N)
  
  # determine N_hat by Bernoulli sampling
  N_hat = stats::rbinom(1, choose(n,r), (N / choose(n,r)))
  
  # Choose randomly N_hat unique subsets with cardinality r of {1,...,n}
  indices = matrix(unlist(random_combs(n,r,N_hat)[[1]]), ncol = r, byrow = TRUE)
  
  # randomly choose minors to test
  ind_minors = random_minors(m,2,nr_minors)
  
  # Calculate matrix H of estimates
  H = H_factors(X, indices, ind_minors)
  p = dim(H)[2]  # total nr of equality constraints

  # Mean and centering
  H_mean = colMeans(H)
  H_centered = t(t(H) - H_mean)
  
  # Compute matrix G
  G = G_factors(X,L=(r-1), ind_minors)
  G_mean = colMeans(G)
  G_centered = t(t(G) - G_mean)
  
  # Diagonal of the approximate variance of H
  cov_H_diag = colSums(H_centered**2) / N_hat
  cov_G_diag = colSums(G_centered**2) / n
  cov_diag = r**2 * cov_G_diag + (n/N) * cov_H_diag
  
  # Vector for standardizing
  standardizer = cov_diag**(-1/2)
  
  # Test statistic
  marginal_stats = abs(sqrt(n) * H_mean)
  test_stat =  max(standardizer * marginal_stats)
  
  # Bootstrap
  bootstrap_res = bootstrap_U(E, r, H_centered, G_centered)
  U_A = bootstrap_res[[1]]
  U_B = bootstrap_res[[2]]
  U = abs(U_A + sqrt(n/N) * U_B)
  U_standardized = t(t(U) * standardizer)
  results = matrixStats::rowMaxs(U_standardized)
  
  # pval
  pval = (1 + sum(results >= test_stat)) / (1+E)
  
  return(list("PVAL"=pval, "TSTAT"=test_stat))
}

