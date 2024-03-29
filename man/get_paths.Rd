% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{get_paths}
\alias{get_paths}
\title{Get the paths between all nodes}
\usage{
get_paths(g)
}
\arguments{
\item{g}{An igraph object that is a tree. It is assumed that the first m nodes correspond to oberseved nodes.}
}
\value{
A list of nested lists. Entry (i,j) contains the path from node i to node j as a list of integers.
}
\description{
Computes the paths between all nodes, observed and unobserved. 
For this task, the function \code{shortest_paths} from the \code{igraph}-package is used.
}
\examples{
vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
tree <- igraph::graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
plot(tree)
res <- get_paths(tree)

# Access path from node 1 to node 2
res[[1]][[2]]
}
