cov_from_graph_large = function(g, paths){
  
  # g is an igraph object with the following attributes
  # V(g)$var: variances for all nodes
  # E(g)$corr: correlation for all edges in (-1,1)
  
  std = sqrt(V(g)$var)
  nr_nodes = length(V(g))
  
  corr = diag(rep(1,nr_nodes))
  for (i in 1:(nr_nodes-1)){
    for (j in (i+1):nr_nodes){
      corr[i,j] = prod(E(g)$corr[paths[[i]][[j]]])
      corr[j,i] = corr[i,j]
    }
  }
  cov = (diag(std) %*% corr %*% diag(std))
  return(cov)
}




update_param = function(g, S){
  
  # g is an igraph object with the following attributes
  # V(g)$var: variances for all nodes
  # E(g)$corr: correlation for all edges in (-1,1)
  # S is large covariance matrix of dimensin equal to total number of nodes
  
  # update correlation of all edges
  edges = get.edgelist(g)
  for (i in 1:dim(edges)[1]){
    from = as.integer(edges[i,][1])
    to = as.integer(edges[i,][2])
    E(g)$corr[i] = S[from, to] / ( sqrt(S[from,from]) * sqrt(S[to,to]) )
    
  }
  
  # update variance of observed nodes
  nr_observed = sum(V(g)$type==1) # always the first ones in vertices
  for (i in 1:nr_observed){
    V(g)$var[i] = S[i,i]
  }
  return(g)
}



# pzwiernik/structuralEM
loglik = function(cov, X){
  
  # cov: SMALL covariance matrix (of observed nodes)
  # X: observed data
  
  n <- dim(X)[1]
  m = dim(X)[2]
  C <- (1/n)*t(X)%*%X
  ll <- - (n/2) * ( m * log(2*pi) + log(det(cov)) + sum(diag(C%*%solve(cov))) )
  
  return(ll)
}



# pzwiernik/structuralEM
E_step <- function(X,S){
  
  # this is the E-step of the algorithm
  # input: the observed data X and the "big" covariance matrix S
  # output: expected sufficient statistics for the model
  
  m <- dim(X)[2] # number of observed nodes
  n <- dim(X)[1] # sample size
  N <- dim(S)[1] # number of nodes in the tree
  
  # in the E-step we compute conditional expectations of sufficient statistics
  # the sufficient statistics are explicitly mentioned in the notes
  # these are the "sample" covariances between H and X (denoted C_XH)
  # and the "sample" variance of H (denoted by C_HH)
  C <- (1/n)*t(X)%*%X
  iSX <- solve(S[1:m,1:m])
  A <- S[(m+1):N,1:m]%*%solve(S[1:m,1:m])%*%C%*%iSX%*%S[1:m,(m+1):N]
  # this is the formula for conditional expectation E( C_HH | X )
  VH_X <- S[(m+1):N,(m+1):N]-S[(m+1):N,1:m]%*%iSX%*%S[1:m,(m+1):N] + A
  V <- matrix(0,nrow=N,ncol=N)
  # this is the formula for conditional expectation E( C_XH | X )
  CXH_X <- C%*%iSX%*%S[1:m,(m+1):N]
  
  # we now create the whole sufficient statistics (including the observed part)
  V[1:m,1:m] <- C
  V[(m+1):N,(m+1):N] <- VH_X
  V[1:m,(m+1):N] <- CXH_X
  V[(m+1):N,1:m] <- t(CXH_X)
  return(V)
}



M_step = function(g, S, paths){
  g = update_param(g,S)
  S = cov_from_graph_large(g, paths)
  return(S)
}



EM = function(X,g,paths,tol=1e-4, maxiter=200){
  m = dim(X)[2] # number of observed nodes
  S = cov_from_graph_large(g, paths)
  l_0 = loglik(S[1:m,1:m], X)
  i=0
  while (i <= maxiter){
    
    # E-step
    S = E_step(X,S)
    
    # M-step
    S = M_step(g,S, paths)
    
    # Evaluate loglik 
    l_1 = loglik(S[1:m,1:m], X)
    print(l_1)
    
    # Check stopping criteria
    if (abs(l_1-l_0) < tol){
      break
    }
    else {
      l_0 = l_1
      i = i+1
    }
  }
  return(list(loglik = l_1, Sigma=S[1:m,1:m]))
}



mle = function(X){
  n = dim(X)[1]
  m = dim(X)[2]
  C = (1/n)*t(X)%*%X
  
  ll = - (n/2) * ( m * log(2*pi) + log(det(C)) + m )
  return(list(loglik=ll, Sigma=C))
}



LR_test = function(X, g, paths){
  
  # returns p value
  LR_statistic = 2*( mle(X)$loglik - EM(X,g,paths)$loglik )
  df = choose((dim(X)[2]+1),2) - ( length(E(g)) + sum(V(g)$type==1) )
  p_value = 1 - stats::pchisq(LR_statistic, df)
  
  return(p_value)
}

