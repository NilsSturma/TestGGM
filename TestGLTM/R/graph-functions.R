findQ = function(g, m){
  
  # g is an igraph object (always first m nodes oberseved)

  sub_sets = CombMSC::subsets(m,4,1:m)
  
  Q = list()
  not_Q = list()
  res = NULL
  
  for (j in 1:nrow(sub_sets)){
    p = sub_sets[j,1]
    others = sub_sets[j,2:4]
    count = 0
    # Find the pair that gives an empty set when the intersection of its component paths is taken
    for (i in 1:3){
      path1 = igraph::shortest_paths(g, from=V(g)[name==p], to=V(g)[name==others[i]], output="epath", weight=NA)
      path2 = igraph::shortest_paths(g, from=V(g)[name==others[-i][1]], to=V(g)[name==others[-i][2]], output="epath", weight=NA)  
      intersection = intersect(unlist(path1$epath), unlist(path2$epath))
      if (length(intersection)==0){
        res = c(p, others[i], others[-i])
        count = count + 1
      }
      if (count>1){
        # If it exists more than one empty intersection of its component paths -> set not in Q
        not_Q = c(not_Q, list(sub_sets[j,1:4]))
        break
      }
      if (i==3){
        Q = c(Q, list(res))
      }
    }
   
  }
  return(list(Q, not_Q))
}




collect_indices <- function(g){
  
  # Determine m (Requirement: First m nodes are always observed.)
  m = sum(V(g)$type==1)
  
  res_findQ = findQ(g, m)
  Q = res_findQ[[1]]
  not_Q = res_findQ[[2]]
  subsets3 = CombMSC::subsets(m,3,1:m)
  
  ind_eq = list()  # Equality constraints (2) (a) and (b)
  ind_ineq1 = list()  # Inequality constraints (1) (a)
  ind_ineq2 = list()  # Inequality constraints (1) (b), (c)
  
  for (i in seq_along(Q)){
    p = Q[[i]][1]
    q = Q[[i]][2]
    r = Q[[i]][3]
    s = Q[[i]][4]
    ind_eq = c(ind_eq, list(c(p,r,q,s,p,s,q,r)))
    ind_ineq2 = c(ind_ineq2, list(c(p,r,q,s,p,q,r,s)))
  }
  for (i in seq_along(not_Q)){
    p = not_Q[[i]][1]
    q = not_Q[[i]][2]
    r = not_Q[[i]][3]
    s = not_Q[[i]][4]
    ind_eq = c(ind_eq, list(c(p,r,q,s,p,s,q,r)), list(c(p,r,q,s,p,q,r,s)))
  }
  for (i in 1:nrow(subsets3)){
    p = subsets3[i,1]
    q = subsets3[i,2]
    r = subsets3[i,3]
    ind_ineq1 = c(ind_ineq1, list(c(p,q,p,r,q,r)))
    ind_ineq2 = c(ind_ineq2, list(c(p,q,q,r,q,q,p,r)), list(c(p,r,q,r,r,r,p,q)), list(c(p,q,p,r,p,p,q,r)))
  }
  
  return(list(ind_eq, ind_ineq1, ind_ineq2))
}




cov_from_graph = function(g){
  
  # g is an igraph object with the following attributes
  # type: 1= observes, 2=unobserved
  # var: variances for all nodes
  # corr: correlation for all edges in (-1,1)
  
  std = sqrt(V(g)$var[V(g)$type==1])
  m = sum(V(g)$type==1)
  
  corr = diag(rep(1,m))
  for (i in 1:m){
    for (j in min((i+1),m):m){
      path = igraph::shortest_paths(g, from=V(g)[name==i], to=V(g)[name==j], output="epath", weight=NA)
      corr[i,j] = prod(E(g)$corr[c(unlist(path$epath))])
      corr[j,i] = corr[i,j]
    }
  }
  cov = (diag(std) %*% corr %*% diag(std))
  return(cov)
}