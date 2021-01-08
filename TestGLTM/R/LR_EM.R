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



# From pzwiernik/structuralEM
loglik = function(cov, X){
  
  # cov: SMALL covariance matrix (of observed nodes)
  # X: observed data
  
  n <- dim(X)[1]
  m = dim(X)[2]
  C <- (1/n)*t(X)%*%X
  ll <- - (n/2) * ( m * log(2*pi) + log(det(cov)) + sum(diag(C%*%solve(cov))) )
  
  return(ll)
}



# From pzwiernik/structuralEM
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


#' Likelihood ratio test for the goodness of fit of a Gaussian latent tree model
#' 
#' Testing the goodness of fit of a given Gaussian latent tree model to observed data.
#' 
#' @param X Matrix with observed data. 
#' Number of columns equal to the number of leaves of the tree (i.e. number of observed variables). 
#' Each row corresponds to one sample.
#' @param g An igraph object that is a tree. It is assumed that the first m nodes correspond to oberseved nodes. 
#' Type 1 indicates that a node is observed. Should be set via \code{V(g)$type==1}.
#' It is assumed that \code{V(g)$var} is the variance of the observed variables and 
#' that \code{E(g)$corr} represents the edge correlations. Should be initialized with starting values.
#' @param paths Nested list with the paths between all nodes. 
#' Should be computed by the function \code{\link{get_paths}}. 
#' This is done outside the LR test to accelerate the computation.
#' @return p-value
#' @examples 
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' plot(tree)
#' 
#' # Set starting values
#' E(tree)$corr = rep(0.7,7)
#' V(tree)$var = rep(1,8)
#' 
#' # Compute all paths
#' paths <- get_paths(tree)
#' 
#' # Call the test
#' LR_test(X, g, paths)
LR_test = function(X, g, paths){
  
  # returns p value
  LR_statistic = 2*( mle(X)$loglik - EM(X,g,paths)$loglik )
  df = choose((dim(X)[2]+1),2) - ( length(E(g)) + sum(V(g)$type==1) )
  p_value = 1 - stats::pchisq(LR_statistic, df)
  
  return(p_value)
}

