library(igraph)
library(CombMSC) # subsets



findQ = function(g){
  
  # g is a graph 
  # V(g)$type indicates whether a node is observed or not (always first m oberseved)
  
  m = sum(V(g)$type==1)
  sub_sets = subsets(m,4,1:m)
  
  Q = list()
  not_Q = list()
  res = NULL
  
  for (j in 1:nrow(sub_sets)){
    p = sub_sets[j,1]
    others = sub_sets[j,2:4]
    count = 0
    # Find the pair that gives an empty set when the intersection of its component paths is taken
    for (i in 1:3){
      path1 = shortest_paths(g, from=V(g)[name==p], to=V(g)[name==others[i]], output="epath", weight=NA)
      path2 = shortest_paths(g, from=V(g)[name==others[-i][1]], to=V(g)[name==others[-i][2]], output="epath", weight=NA)  
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
      path = shortest_paths(g, from=V(g)[name==i], to=V(g)[name==j], output="epath", weight=NA)
      corr[i,j] = prod(E(g)$corr[c(unlist(path$epath))])
      corr[j,i] = corr[i,j]
    }
  }
  cov = (diag(std) %*% corr %*% diag(std))
  return(cov)
}