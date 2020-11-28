library(StructuralEM)
#library(igraph)
library(ape)

m <- 5 # number of leaves
n <- 30 # sample size


Ttr <- rtree(m,rooted=FALSE,br=rep(-log(.6),7))
Ttr$tip.label <- 1:m # label leaves by 1:m
plot(Ttr)

Str <- get.corr0(Ttr) # correlation matrix all variable
Rtr <- Str[1:m,1:m]  # correlation matrix observed variables


dat <- mvrnorm(n, rep(0,m), Rtr) # randomly sample observed variables


D <- get.dist(cor(dat))  # estimated distance matrix (from estmated correlation matrix)

T0 <- fastme.bal(D)
T0$edge.length <- (T0$edge.length>0)*T0$edge.length
res <- strEM(dat,T0,tol=1e-6)
