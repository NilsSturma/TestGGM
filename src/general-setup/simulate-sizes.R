library(foreach)
library(doParallel)
library(MASS) #mvrnorm
source("graph-functions.R")
source("tests.R")


#################
# Create a tree #
#################

# !!! Important: observed variables = leaves !!!


# Quinted tree
colors <- c("tomato", "gray50")
vertices <- data.frame(name=seq(1,8),
                     type=c(rep(1,5), rep(2,3)), #1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                     color=colors[c(1,1,1,1,1,2,2,2)])
edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
plot(g)
V(g)$var = rep(2,8)
E(g)$corr = rep(0.5,7)
cov = cov_from_graph(g)

# Star tree
m = 10
colors <- c("tomato", "gray50")
vertices <- data.frame(name=seq(1,m+1), type=c(rep(1,m), 2), color=colors[c(rep(1,m), 2)])
edges <- data.frame(from=seq(1,m), to=rep((m+1),m))
g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
plot(g)
# Setup1
V(g)$var = rep(2,11)
E(g)$corr = rep(sqrt(0.5),10)
cov = cov_from_graph(g)
# # Setup2 (Sample in loop !!!)
# beta = c(10,10, rnorm((m-2),0,0.2))
# V(g)$var = c(beta**2+rep(1/3,m),1)
# E(g)$corr = beta / sqrt(beta**2+rep(1/3,m)) 
# cov = cov_from_graph(g)
# cov


#########################
# Determine Q and not_Q #
#########################

res = findQ(g)
Q = res[[1]]
not_Q = res[[2]]


###################
# Run experiments #
###################

# Set variables
n = 500
B = 3  # just for "1-dependent"
E = 1000
alphas = seq(0.01, 0.99, 0.01)
nr_exp = 500
method="1-dependent"

save=TRUE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "quinted-tree_", method, "_n=", n, sep="")


# Emprirical size for each alpha
cores = detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)
sizes <- foreach(nr = 1:nr_exp, .combine=rbind, .packages=c("MASS")) %dopar% {
  # Generate n independent data sets
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  # Call the test
  result = test_equality_constraints(X, Q, not_Q, method=method, B=B, E=E, alphas=alphas)
  as.numeric(result)
}
sizes = colMeans(sizes)
stopCluster(cl)



#########################
# Plot and save results #
#########################

subtitle = paste("Quinted tree - n=", n, sep="")
title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments", sep="")

# Plot
if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste(name, ".pdf", sep="")
  name_rds = paste(name, ".rds", sep="")
  saveRDS(sizes, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}

plot(alphas, sizes, 
     xlab="Nominal level", ylab="Emprical test size", main=title, #sub=subtitle,
     type="p", pch=1)
abline(coef = c(0,1))

if (save){
  dev.off() # close pdf file
}



















##########
# Others #
##########

# # Random multivariate normal distribution
# A = matrix(rnorm(25,0,1), nrow=5)
# Sigma = A %*% t(A)
# X = mvrnorm(500, mu=rep(0,5), Sigma=Sigma)
