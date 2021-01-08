#' Samples data from a Gaussian latent tree model
#' 
#' @param tree An igraph object that is a tree. 
#' It is assumed that the first \code{m} nodes correspond to oberseved nodes and 
#' that model parameters are given in the following way:
#' \code{V(g)$var} is the variance of the observed variables and 
#' that \code{E(g)$corr} are the edge correlations. 
#' @param m Integer, number of observed nodes.
#' @param n Integer, number of samples being sampled.
#' @return Matrix with sampled data with m columns and n rows. 
#' Each row corresponds to one sample.
#' @examples
#' # Create tree
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' 
#' # Sample data from tree
#' V(tree)$var = rep(1,8)
#' E(tree)$corr = rep(0.7,7)
#' sample_from_tree(tree, m=5, n=500)
sample_from_tree <- function(tree, m, n){
  
  # TODO: Add checks for correct inputs
  # TODO: Add error handling
  
  paths = get_paths(tree)
  cov = cov_from_graph(tree, m, paths)
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  return(X)
}





#' Tests the goodness-of-fit of a Gaussian latent tree model
#' 
#' This function tests the goodness-of-fit of a given Gaussian latent tree model to observed data.
#' It supports all Gaussian latent tree models where observed variables are restricted to 
#' be the leaves of the tree. Four different test strategies are implmented. 
#' One is the likelihood ratio test, the other three are algebraic tests. 
#' In the latter case the test statistic is formed as the maximum of unbiased estimates of the 
#' polynomials determining the semialgebraic set that defines the parameter space of the model. 
#' A Gaussian multiplier bootstrap is used to estimate the limiting distribution of the 
#' test statistic and to compute the p-value of the test.
#' 
#' @param X Matrix with observed data. Number of columns has to be equal to the number of 
#' leaves of the tree (i.e. number of observed variables). Each row corresponds to one sample.
#' @param tree An igraph object that is a tree. It is assumed that the first \code{m} nodes correspond to oberseved nodes.
#' @param m Integer, number of observed nodes.
#' @param test_strategy String, determines the test that is applied. Has to be one 
#' out of \code{c("LR", "grouping", "run-over", "U-stat")}. Default is \code{"gouping"}.
#' @param E Integer, number of bootstrap iterations. This has no effect if \code{test_strategy="LR"}.
#' @param B Integer, batch size for the estimate of the covariance matrix. 
#' Only relevant for the \code{"run-over"} test.
#' @param N Integer, computational budget parameter. Only relevant for the \code{"U-stat"} test.
#' @param only_equalities Logical. Should the test incorporate only equality constraints? 
#' Default \code{FALSE}. This has no effect if \code{test_strategy="LR"}.
#' @param nr_4 Number of considered subsets of size 4. This is optional. 
#' Useful for very high-dimensional models. Default \code{NULL} hence all subsets are considered. 
#' If a number is given, subsets are choosen randomly and only the corresponding constraints are considered.
#' @param nr_3 Number of considered subsets of size 3. This is optional. 
#' Useful for very high-dimensional models. Default \code{NULL} hence all subsets are considered. 
#' If a number is given, subsets are choosen randomly and only the corresponding constraints are considered.
#' 
#' @return Named list with two entries: Test statistic (\code{TSTAT}) and p-value (\code{PVAL}).
#' 
#' @examples
#' # Create tree
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' 
#' # Sample data from tree
#' V(tree)$var = rep(1,8)
#' E(tree)$corr = rep(0.7,7)
#' X = sample_from_tree(tree, m=5, n=500)
#' 
#' # Goodness of fit test
#' gltmtest(X, tree, m=5, test_strategy="grouping")
gltmtest <- function(X, tree, m, 
                    test_strategy="grouping", 
                    E=1000, 
                    B=5, 
                    N=5000, 
                    only_equalities=FALSE, 
                    nr_4 = NULL, 
                    nr_3=NULL){
  
  # TODO: Add checks for correct inputs
  # TODO: Add error handling
  
  # Collect all paths in the given tree
  paths = get_paths(tree)
  
  # Determine polynomial constraints to be tested
  if (test_strategy!="LR"){
    res = collect_indices(tree, m, nr_4, nr_3)
    ind_eq = res$ind_eq
    ind_ineq1 = res$ind_ineq1
    ind_ineq2 = res$ind_ineq2
    p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
    if (only_equalities){
      ind_ineq1 = NULL
      ind_ineq2 = NULL
    }
  }
  
  # Call the test
  if (test_strategy=="LR"){
    res = LR_test(X,tree,paths)
  } else if (test_strategy=="grouping"){
    res = test_grouping(X, ind_eq, ind_ineq1, ind_ineq2, E=E)
  } else if (test_strategy=="run-over"){
    res = test_run_over(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E)
  } else if (test_strategy=="U-stat"){
    res = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E)
  } 
  return(res)
}