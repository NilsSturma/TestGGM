library(foreach)
library(doParallel)
library(MASS) #mvrnorm
library(igraph)
library(TestGLTM)


source("simulations/utils.R") # TODO: add these functions to package

#################
# Set variables #
#################

# General
n_range = c(100, 250, 500, 1000)
#n = 1000
E = 1000
nr_exp = 500
alphas = seq(0.01, 0.99, 0.01)

# Test strategy
test_strategy="U-stat-deg"  # "grouping", "run-over", "U-stat", "LR", "U-stat-deg"
B = 5  # just for test_strategy=="run-over" (5 works best for setup 1 after doing some experiments)
N = 5000 # the more, the better but not feasible
#N_range = c(2*n, 5*n, round(n**1.5), round(n**1.8), round(n**2))

# Tree
tree = "cat_binary"  # "star_tree", "cat_binary"
m = 20  # (star_tree)
setup = 1  # (star_tree)


# Saving
save=TRUE


###################################
# Create tree and collect indices #
###################################

if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="cat_binary"){
  g = cat_binary()
} 

plot(g)

res = collect_indices(g)
ind_eq = matrix(unlist(res[[1]]), ncol = 8, byrow = TRUE)
ind_ineq1 = matrix(unlist(res[[2]]), ncol = 6, byrow = TRUE)
ind_ineq2 = matrix(unlist(res[[3]]), ncol = 8, byrow = TRUE)
p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
print(p)


#######################################
# Test wrapper function (single test) #
#######################################
# n = 1000
# X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
# bootstrap_test(X, g, alpha=0.05, method="run-over", B=5)


for (n in n_range){
  print(paste("n=",n ,sep=""))
  ###############################################
  # Compute empirical test sizes for all alphas #
  ###############################################
  
  cores = 20 # detectCores()
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  results <- foreach(nr = 1:nr_exp, 
                     .combine=rbind, 
                     .errorhandling="remove",   # "pass", "stop", "remove"
                     .packages=c("MASS", "TestGLTM", "igraph", "stats")) %dopar% {
    
    if((nr%%10) == 0){
      print(nr)
    }
    warnings()
    
    # Generate n independent data sets depending on setup
    if (tree=="star_tree"){
      cov = cov_from_star_tree(g, setup=setup, m=m)
    } else if (tree=="cat_binary"){
      V(g)$var = rep(1,38)
      E(g)$corr = rep(0.7,37)
      cov = cov_from_graph(g)
    }
    
    X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
    
    # Call the test
    if (test_strategy=="LR"){
      if (tree=="star_tree"){
        res = factanal(X, 1)
        result = res[["PVAL"]] <= alphas # result: TRUE = rejected
      } else if (tree=="cat_binary"){
        result = LR_test(X,g) <= alphas # result: TRUE = rejected
      }
    } else if (test_strategy=="grouping"){
      result = test_grouping(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alphas)
    } else if (test_strategy=="run-over"){
      result = test_run_over(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E, alphas=alphas)
    } else if (test_strategy=="U-stat"){
      result = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E, alphas=alphas)
    } else if (test_strategy=="U-stat-deg"){
      result = test_U_stat_degenerate(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E, alphas=alphas)
    }
    
    
    result = as.numeric(result)
  }
  
  # Check if every task was successful
  if (dim(results)[1] != nr_exp){
    print(paste("ERROR - ", (nr_exp-dim(results)[1]), " tasks where not succesful.", sep=""))
  }
  
  sizes = colMeans(results)
  stopCluster(cl)
  
  
  
  #########################
  # Plot and save results #
  #########################
  
  if (tree=="star_tree"){
    name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "star-tree_setup=", setup, "_n=", n, "_m=", m, sep="")
    subtitle = paste("Star tree n=", n, " m=", m," strategy=", test_strategy, sep="")
    title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments. \n Star tree - setup ", setup, sep="")
  } else if (tree=="cat_binary"){
    name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "caterpillar", "_n=", n,  sep="")
    subtitle = paste("Caterpillar tree n=", n, " m=20 strategy=", test_strategy, sep="")
    title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments. \n Caterpillar tree", sep="")
  }
  
  #name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "star-tree_setup=", setup, "_N=", N, sep="")
  
  
  
  # Plot
  if (save){
    # use "./img/name.png" to save in subdirectory
    name_pdf = paste("./results/", tree, "/", test_strategy, "/sizes/", name, ".pdf", sep="")
    name_rds = paste("./results/", tree, "/", test_strategy, "/sizes/", name, ".rds", sep="")
    saveRDS(sizes, file = name_rds) # read with readRDS()
    pdf(name_pdf) # create pdf file
  }
  
  plot(alphas, sizes, 
       xlab="Nominal level", ylab="Emprical test size", main=title, sub=subtitle,
       type="p", pch=1)
  abline(coef = c(0,1))
  if (save){
    dev.off() # close pdf file
  }
}
