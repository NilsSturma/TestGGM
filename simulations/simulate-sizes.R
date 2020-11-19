library(foreach)
library(doParallel)
library(MASS) #mvrnorm
library(igraph)
library(TestGLTM)

setwd("src/general-setup")

#################
# Create a tree #
#################

# TODO: put the following three functions in the package

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


# Create an igraph object that is star tree
star_tree <- function(m=10){
  colors <- c("tomato", "gray50")
  vertices <- data.frame(name=seq(1,m+1), 
                         type=c(rep(1,m), 2), # 1=observed, 2=latent -> always first m nodes should be observed (=leaves)
                         color=colors[c(rep(1,m), 2)])
  edges <- data.frame(from=seq(1,m), to=rep((m+1),m))
  g <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
  return(g)
}

# Sample covariance matrix from star tree depending on setup
cov_from_star_tree <- function(g, setup=1){
  if (setup==1){
    V(g)$var = rep(2,11)
    # for unobserved node: variance does not matter
    E(g)$corr = rep(sqrt(0.5),10)
    cov = cov_from_graph(g)
  } else if (setup==2){
    beta = c(10,10, rnorm((m-2),0,0.2))
    V(g)$var = c(beta**2+rep(1/3,m),1)
    E(g)$corr = beta / sqrt(beta**2+rep(1/3,m))
    cov = cov_from_graph(g)
  }
  return(cov)
}








#################
# Set variables #
#################

# Set variables
n = 100000
m = 10
E = 1000
setup = 2
alphas = seq(0.01, 0.99, 0.01)
nr_exp = 500

test_strategy="run-over"  # "run-over", "grouping"
B = 5  # just for test_strategy=="run-over" (5 works best for setup 1 after doing some experiments)

save=FALSE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "quinted-tree_", method, "_n=", n, sep="")

################
# Create graph #
################
g = star_tree(m)
plot(g)
# g = quinted_tree()
# plot(g)
# V(g)$var = rep(2,8)
# E(g)$corr = rep(0.5,7)
# cov = cov_from_graph(g)


# Collect indices
res = collect_indices(g)
ind_eq = matrix(unlist(res[[1]]), ncol = 8, byrow = TRUE)
ind_ineq1 = matrix(unlist(res[[2]]), ncol = 6, byrow = TRUE)
ind_ineq2 = matrix(unlist(res[[3]]), ncol = 8, byrow = TRUE)





#######################################
# Test wrapper function (single test) #
#######################################
# n = 1000
# X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
# bootstrap_test(X, g, alpha=0.05, method="run-over", B=5)



###############################################
# Compute empirical test sizes for all alphas #
###############################################

cores = 20  # detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

results <- foreach(nr = 1:nr_exp, .combine=rbind, .packages=c("MASS", "TestGLTM", "igraph")) %dopar% {
  
  if((nr%%10) == 0){
    print(nr)
  }
  
  # Generate n independent data sets depending on setup
  cov = cov_from_star_tree(g, setup=setup)
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  
  # Call the test
  if (test_strategy=="grouping"){
    result = test_independent(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alphas)
  } else if (test_strategy=="run-over"){
    result = test_m_dep(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E, alphas=alphas)
  }
  result = as.numeric(result)
}

sizes = colMeans(results)
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
