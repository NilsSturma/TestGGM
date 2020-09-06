library(MASS) #mvrnorm
source("graph-functions.R")
source("general-tests.R")

#################
# Create a tree #
#################

# Important: observed variables = leaves
colors <- c("tomato", "gray50")
vertices <- data.frame(name=seq(1,8),
                     type=c(rep(1,5), rep(2,3)), #1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                     color=colors[c(1,1,1,1,1,2,2,2)]) 
edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
plot(g)



# Find set Q for this tree
res = findQ(g)
Q = res[[1]]
not_Q = res[[2]]


##########################
# Sample from given tree #
##########################

V(g)$var = rep(2,8)
E(g)$corr = rep(0.5,7)
cov = cov_from_graph(g)
X = mvrnorm(500, mu=rep(0,5), Sigma=cov)

#########################
# Determine Q and not_Q #
#########################

res = findQ(g)
Q = res[[1]]
not_Q = res[[2]]

###################
# Run experiments #
###################
bla = test_equality_constraints(X, Q, not_Q, alphas=0.5)
bla







##########
# Others #
##########

# # Random multivariate normal distribution
# A = matrix(rnorm(25,0,1), nrow=5)
# Sigma = A %*% t(A)
# X = mvrnorm(500, mu=rep(0,5), Sigma=Sigma)
