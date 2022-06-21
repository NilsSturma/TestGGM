# Create an igraph object that is the quintet tree
quinted_tree <- function(){
  colors <- c("tomato", "gray50")
  vertices <- data.frame(name=seq(1,8),
                         type=c(rep(1,5), rep(2,3)), # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[c(1,1,1,1,1,2,2,2)])
  edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
  g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}



# Create an igraph object that is a binary rooted tree with 12 observed nodes (i.e 10 hidden nodes, 21 edges)
binary_rooted <- function(){
  colors <- c("tomato", "gray50")
  name = seq(1,22)
  type = c(rep(1,12), rep(2,10))
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=c(22,22,22,19,19,20,20,21,21,13,13,14,14,15,15,16,16,17,17,18,18), 
                      to=c(19,20,21,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12))
  g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}


# Create an igraph object that is a binary caterpillar tree with m observed variables
cat_binary <- function(m=10){
  colors <- c("tomato", "gray50")
  name = seq(1,(m+(m-2)))
  type = c(rep(1,m), rep(2,(m-2)))
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=c(seq((m+1),(m+(m-2))),(m+1),(m+(m-2)),seq((m+1),(m+(m-2)-1))),
                      to=c(seq(1,(m-2)),(m-1),(m),seq((m+2),(m+(m-2)))) )
  g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}



# Create an igraph object that is star tree with m observed variables
star_tree <- function(m=10){
  colors <- c("tomato", "gray50")
  name = seq(1,m+1)
  type = c(rep(1,m), 2)
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=seq(1,m), to=rep((m+1),m))
  g <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}



# Sample covariance matrix from star tree depending on setup
cov_from_star_tree <- function(g, paths, setup=1, m=10){
  if (setup==1){
    igraph::V(g)$var = c(rep(2,m),1)
    igraph::E(g)$corr = rep(sqrt(0.5),m)
    cov = cov_from_graph(g, m, paths)
  } else if (setup==2){
    beta = c(10,10, stats::rnorm((m-2),0,0.2))
    igraph::V(g)$var = c(beta**2+rep(1/3,m),1)
    igraph::E(g)$corr = beta / sqrt(beta**2+rep(1/3,m))
    cov = cov_from_graph(g, m, paths)
  } else if (setup==3){
    beta = stats::rnorm(m,0,0.2)
    igraph::V(g)$var = c(beta**2+rep(1/3,m),1)
    igraph::E(g)$corr = beta / sqrt(beta**2+rep(1/3,m))
    cov = cov_from_graph(g, m, paths)
  }
  return(cov)
}

