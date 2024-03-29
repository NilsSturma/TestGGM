#' Determine the set Q
#' 
#' For all subsets of size 4 of the observed nodes of a tree, the function determines if a subset is an element of Q or not.
#' 
#' @param g An igraph object that is a tree. It is assumed that the first m nodes correspond to observed nodes.
#' @param m Integer, number of observed nodes.
#' @param nr Integer, number of subsets that are considered. Note that this is optional. If \code{NULL}, all subsets are considered.
#' @return A list with two entries is returned. Both entries are matrices with 4 columns. 
#' The first matrix are all subsets that belong to Q, the second are all subsets that do not belong to Q. 
#' Each row represents one subset. For the matrix Q the order of the elements in a row is important. 
#' A row with entries p,q,r,s means that the element {p,q}|{r,s} is in Q.
#' @examples
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' plot(tree)
#' findQ(tree, 5)
#' @export
findQ = function(g, m, nr=NULL){
  
  if (is.null(nr)){
    sub_sets = t(utils::combn(m,4))
  } else {
    sub_sets = matrix(unlist(random_combs(m,4,nr)[[1]]), ncol = 4, byrow = TRUE)
  }
  
  Q = list()
  not_Q = list()
  res = NULL
  
  for (j in 1:nrow(sub_sets)){
    p = sub_sets[j,1]
    others = sub_sets[j,2:4]
    count = 0
    # Find the pair that gives an empty set when the intersection of its component paths is taken
    for (i in 1:3){
      path1 = igraph::shortest_paths(g, from=p, to=others[i], output="epath", weight=NA)
      path2 = igraph::shortest_paths(g, from=others[-i][1], to=others[-i][2], output="epath", weight=NA)  
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







#' Collect polynomial constraints that characterize the semialgebraic variety of a given tree
#' 
#' The polynomials are represented as a list of indices. Example:
#' If one polynomial is \eqn{\sigma_{pq} \sigma_{rs} - \sigma_{pr} \sigma_{qs}} then the indices will be \code{p,q,r,s,p,r,q,s}. 
#' 
#' @param g An igraph object that is a tree. It is assumed that the first m nodes correspond to observed nodes.
#' @param m Integer, number of observed nodes.
#' @param nr_4 Number of considered subsets of size 4. This is optional. 
#' If \code{NULL}, all subsets are considered. 
#' If a number is given, subsets are chosen randomly.
#' @param nr_3 Number of considered subsets of size 3. This is optional. 
#' If \code{NULL}, all subsets are considered. 
#' If a number is given, subsets are chosen randomly.
#' @return A list with three entries is returned. All entries are matrices. 
#' The first matrix contains the equality constraints and has 8 columns. 
#' The second matrix contains all inequality constraints where only 6 indices are necessary (i.e. it has 6 columns). 
#' The third matrix contains the inequality constraints where 8 indices are necessary (i.e. it has 8 columns).
#' @examples
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' plot(tree)
#' collect_indices(tree, 5)
#' @export
collect_indices <- function(g, m, nr_4=NULL, nr_3=NULL){
  
  res_findQ = findQ(g, m, nr_4)
  Q = res_findQ[[1]]
  not_Q = res_findQ[[2]]
  
  if (is.null(nr_3)){
    sub_sets3 = t(utils::combn(m,3))
  } else {
    sub_sets3 = matrix(unlist(random_combs(m,3,nr_3)[[1]]), ncol = 3, byrow = TRUE)
  }
  
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
  for (i in 1:nrow(sub_sets3)){
    p = sub_sets3[i,1]
    q = sub_sets3[i,2]
    r = sub_sets3[i,3]
    ind_ineq1 = c(ind_ineq1, list(c(p,q,p,r,q,r)))
    ind_ineq2 = c(ind_ineq2, list(c(p,q,q,r,q,q,p,r)), list(c(p,r,q,r,r,r,p,q)), list(c(p,q,p,r,p,p,q,r)))
  }
  
  return(list("ind_eq" = matrix(unlist(ind_eq), ncol = 8, byrow = TRUE),
              "ind_ineq1" = matrix(unlist(ind_ineq1), ncol = 6, byrow = TRUE), 
              "ind_ineq2" = matrix(unlist(ind_ineq2), ncol = 8, byrow = TRUE)))
}



#' Get the paths between all nodes
#'
#' Computes the paths between all nodes, observed and unobserved. 
#' For this task, the function \code{shortest_paths} from the \code{igraph}-package is used.
#'
#' @param g An igraph object that is a tree. It is assumed that the first m nodes correspond to oberseved nodes.
#' @return A list of nested lists. Entry (i,j) contains the path from node i to node j as a list of integers.
#' 
#' @examples
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' plot(tree)
#' res <- get_paths(tree)
#' 
#' # Access path from node 1 to node 2
#' res[[1]][[2]]
#' @export
get_paths = function(g){
  nr_nodes = length(igraph::V(g))
  paths = list()
  for (i in 1:nr_nodes){
    p_list = list()
    for (j in 1:nr_nodes){
      p = igraph::shortest_paths(g, from=i, to=j, output="epath", weight=NA)
      p = c(unlist(p$epath))
      p_list = c(p_list, list(p))
    }
    paths[[i]] = p_list
  }
  return(paths)
}






#' Covariance matrix of observed variables 
#' 
#' Given a tree and parameters (variance and edge correlations) the covariance matrix of the corresponding multivariate Gaussian distribution is computed.
#'
#' @param g An igraph object that is a tree. It is assumed that the first m nodes correspond to oberseved nodes. 
#' Type 1 indicates that a node is observed. Should be set via \code{V(g)$type==1}.
#' It is assumed that \code{V(g)$var} is the variance of the observed variables and that \code{E(g)$corr} represents the edge correlations.
#' @param m Integer, number of observed nodes.
#' @param paths Nested list with the paths between all nodes. 
#' Should be computed with the function \code{\link{get_paths}}. 
#' @return Covarianz matrix.
#' @examples 
#' vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
#' edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
#' tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
#' plot(tree)
#' 
#' # Set parameters
#' igraph::E(tree)$corr = rep(0.7,7)
#' igraph::V(tree)$var = rep(1,8)
#' 
#' # Compute all paths
#' paths <- get_paths(tree)
#' 
#' # Call function
#' cov_from_graph(tree, 5, paths)
#' @export
cov_from_graph = function(g, m, paths){
  
  std = sqrt(igraph::V(g)$var[igraph::V(g)$type==1])
  m = sum(igraph::V(g)$type==1)
  
  corr = diag(rep(1,m))
  for (i in 1:(m-1)){
    for (j in (i+1):m){
      corr[i,j] = prod(igraph::E(g)$corr[paths[[i]][[j]]])
      corr[j,i] = corr[i,j]
    }
  }
  cov = (diag(std) %*% corr %*% diag(std))
  return(cov)
}


