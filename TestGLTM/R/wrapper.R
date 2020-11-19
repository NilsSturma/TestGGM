bootstrap_test <- function(X, tree, method="grouping", E=1000, alpha=0.05, B=5){
  # TODO: Add checks
  # TODO: Change this. Should return value of test statistic, estimated quantile, etc.
  # Returns: TRUE if nullhypothesis (data X comes from given tree) is rejected
  
  # Collect indices
  indices = collect_indices(tree)
  ind_eq = matrix(unlist(indices[[1]]), ncol = 8, byrow = TRUE)
  ind_ineq1 = matrix(unlist(indices[[2]]), ncol = 6, byrow = TRUE)
  ind_ineq2 = matrix(unlist(indices[[3]]), ncol = 8, byrow = TRUE)
  
  # Call test
  if (method=="grouping"){
    is_rejected = test_independent(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alpha)
  } else if (method=="run-over"){
    is_rejected = test_m_dep(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=1000, alphas=alpha)
  } 
  return(is_rejected)
}