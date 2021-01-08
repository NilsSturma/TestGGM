gofgltm <- function(X, tree, m, 
                    method="grouping", 
                    E=1000, 
                    B=5, 
                    N=5000, 
                    only_equalities=FALSE, 
                    nr_4 = NULL, 
                    nr_3=NULL){
  
  # TODO: Add checks

  paths = get_paths(tree)
  
  res = collect_indices(tree, m, nr_4, nr_3)
  ind_eq = res$ind_eq
  ind_ineq1 = res$ind_ineq1
  ind_ineq2 = res$ind_ineq2
  p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
  if (only_equalities){
    ind_ineq1 = NULL
    ind_ineq2 = NULL
    p = dim(ind_eq)[1]
  }
  
  # Call the test
  if (test_strategy=="LR"){
    if (tree=="star_tree"){
      res = factanal(X, 1)
    } else if (tree=="cat_binary"){
      res = LR_test(X,tree,paths)
    }
  } else if (test_strategy=="grouping"){
    res = test_grouping(X, ind_eq, ind_ineq1, ind_ineq2, E=E)
  } else if (test_strategy=="run-over"){
    res = test_run_over(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E)
  } else if (test_strategy=="U-stat"){
    res = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E)
  } 
  return(res)
}