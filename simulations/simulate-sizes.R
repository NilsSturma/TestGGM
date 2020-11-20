library(foreach)
library(doParallel)
library(MASS) #mvrnorm
library(igraph)
library(TestGLTM)

setwd("../simulations")
source("utils.R")

#################
# Set variables #
#################

# General
n = 1000
E = 1000
nr_exp = 500
alphas = seq(0.01, 0.99, 0.01)


# Test strategy
test_strategy="symmetric"  # "run-over", "grouping", "two-step", "symmetric"
B = 5  # just for test_strategy=="run-over" (5 works best for setup 1 after doing some experiments)
beta = 0.001  # just for test_strategy=="two-step"

# Tree
tree = "star_tree"  # "star_tree", "quinted_tree", "binary_rooted"
m = 10  # (star_tree)
setup = 2  # (star_tree)



# Saving
save=FALSE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "quinted-tree_", method, "_n=", n, sep="")



###################################
# Create tree and collect indices #
###################################

if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="binary_rooted"){
  g = binary_rooted()
} else if (tree=="quinted_tree"){
  g = quinted_tree()
}

plot(g)

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
  warnings()
  
  # Generate n independent data sets depending on setup
  if (tree=="star_tree"){
    cov = cov_from_star_tree(g, setup=setup, m=m)
  } else if (tree=="binary_rooted"){
    V(g)$var = rep(2,22)
    E(g)$corr = rep(0.5,21)
    cov = cov_from_graph(g)
  } else if (tree=="quinted_tree"){
    V(g)$var = rep(2,8)
    E(g)$corr = rep(0.5,7)
    cov = cov_from_graph(g)
  }
 
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  
  # Call the test
  if (test_strategy=="grouping"){
    result = test_independent(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alphas)
  } else if (test_strategy=="run-over"){
    result = test_m_dep(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E, alphas=alphas)
  } else if (test_strategy=="two-step"){
    result = test_two_step(X, ind_eq, ind_ineq1, ind_ineq2, E=E, beta=beta, alphas=alphas)
  } else if (test_strategy=="symmetric"){
    result = test_symmetric(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alphas)
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
