library(TestGLTM)
library(igraph)


# Create an igraph object that is the quinted tree
quinted_tree <- function(){
  colors <- c("tomato", "gray50")
  vertices <- data.frame(name=seq(1,8),
                         type=c(rep(1,5), rep(2,3)), # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[c(1,1,1,1,1,2,2,2)])
  edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}



# Create an igraph object that is a binar rooted tree with 12 observed nodes (i.e 10 hidden nodes, 21 edges)
binary_rooted <- function(){
  colors <- c("tomato", "gray50")
  name = seq(1,22)
  type = c(rep(1,12), rep(2,10))
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=c(22,22,22,19,19,20,20,21,21,13,13,14,14,15,15,16,16,17,17,18,18), 
                      to=c(19,20,21,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}

cat_binary <- function(){
  colors <- c("tomato", "gray50")
  name = seq(1,38)
  type = c(rep(1,20), rep(2,18))
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=c(seq(21,38),21,38,seq(21,37)),
                      to=c(seq(1,18),19,20,seq(22,38)))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}

cat <- function(){
  colors <- c("tomato", "gray50")
  name = seq(1,29)
  type = c(rep(1,20), rep(2,9))
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=c(1,2,seq(3,11),seq(12,20),seq(21,28)),
                      to=c(21,29,seq(21,29),seq(21,29),seq(22,29)))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}




# Create an igraph object that is star tree
star_tree <- function(m=10){
  colors <- c("tomato", "gray50")
  name = seq(1,m+1)
  type = c(rep(1,m), 2)
  vertices <- data.frame(name=name, 
                         type=type, # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[type])
  edges <- data.frame(from=seq(1,m), to=rep((m+1),m))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}



# Sample covariance matrix from star tree depending on setup
cov_from_star_tree <- function(g, setup=1, m=10){
  if (setup==1){
    V(g)$var = rep(2,m+1)
    # for unobserved node: variance does not matter
    E(g)$corr = rep(sqrt(0.5),m)
    cov = cov_from_graph(g)
  } else if (setup==2){
    beta = c(10,10, rnorm((m-2),0,0.2))
    V(g)$var = c(beta**2+rep(1/3,m),1)
    E(g)$corr = beta / sqrt(beta**2+rep(1/3,m))
    cov = cov_from_graph(g)
  }
  return(cov)
}

