library(foreach)
library(doParallel)
library(MASS)
library(igraph)
library(TestGLTM)

setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")

#################
# Set variables #
#################

# General
n_range = c(500)
E = 1000
nr_exp = 500
alphas = seq(0.01, 0.99, 0.01)

# Test strategy
strategies = c("LR")  # Possible: "grouping", "run-over", "U-stat", "LR"
B = 5  # only relevant if test_strategy=="run-over" 
N = 5000  # only relevant if test_strategy=="U-stat"

# Tree
tree = "star_tree"  # Possible: "star_tree", "cat_binary"
m = 20
setup = 1  # only relevant if tree=="star_tree"

# High dimensionality?
nr_4 = NULL  # 5000, NULL
nr_3 = NULL  # 125, NULL

# Test only equalities?
only_equalities = FALSE

# Saving
save=TRUE



###################################
# Create tree and collect indices #
###################################

# Create tree
if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="cat_binary"){
  g = cat_binary(m)
} 

# Plot tree
plot(g)

# Save all paths between all nodes in the tree (doing this just once reduces computational time)
paths = get_paths(g)

# Collect the representations of the polynomials that have to be tested
res = collect_indices(g, m, nr_4, nr_3)
ind_eq = res$ind_eq
ind_ineq1 = res$ind_ineq1
ind_ineq2 = res$ind_ineq2
p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
if (only_equalities){
  ind_ineq1 = NULL
  ind_ineq2 = NULL
  p = dim(ind_eq)[1]
}

# Check the dimension
print(p)



#################################
# Simulate empirical test sizes #
#################################

for (test_strategy in strategies){
  for (n in n_range){
    
    print(paste("strategy=",test_strategy ,sep=""))
    print(paste("n=",n ,sep=""))
    
    # Initialize the cluster
    cores = 20 #detectCores()
    cl <- makeCluster(cores, outfile = "")
    registerDoParallel(cl)
    
    # Main loop to compute test sizes
    results <- foreach(nr = 1:nr_exp, 
                       .combine=rbind, 
                       .errorhandling="remove",   # "pass", "stop", "remove"
                       .packages=c("MASS", "TestGLTM", "igraph", "stats")) %dopar% {
      
      # Print some info
      if((nr%%10) == 0){print(nr)}
      warnings()
      
      # Generate n independent data sets depending on setup
      if (tree=="star_tree"){
        cov = cov_from_star_tree(g, paths, setup=setup, m=m)
      } else if (tree=="cat_binary"){
        V(g)$var = rep(1,(m+(m-2)))
        E(g)$corr = rep(0.7,(m+(m-3)))
        cov = cov_from_graph(g, m, paths)
      }
      X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
      
      # Call the test
      if (test_strategy=="LR"){
        if (tree=="star_tree"){
          res = factanal(X, 1)
        } else if (tree=="cat_binary"){
          res = LR_test(X,g,paths)
        }
      } else if (test_strategy=="grouping"){
        res = test_grouping(X, ind_eq, ind_ineq1, ind_ineq2, E=E)
      } else if (test_strategy=="run-over"){
        res = test_run_over(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E)
      } else if (test_strategy=="U-stat"){
        res = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E)
      } 
      
      # Rejected?
      result = res$PVAL <= alphas # result: TRUE = rejected
      result = as.numeric(result)
    }
    
    # Check if every task was successful
    if (dim(results)[1] != nr_exp){
      print(paste("ERROR - ", (nr_exp-dim(results)[1]), " tasks where not succesful.", sep=""))
    }
    
    # Compute the emprical sizes
    sizes = colMeans(results)
    
    # Stop the cluster
    stopCluster(cl)
    
    # Plot and save results
    # Create meta information for the plot and the files to be saved
    if (only_equalities){
      prefix = "only_equalities"
    } else {
      prefix = format(Sys.time(), "%Y-%m-%d-%H-%M")
    }
    if (tree=="star_tree"){
      name = paste(prefix, "_", "star-tree_setup=", setup, "_n=", n, "_m=", m, sep="")
      subtitle = paste("Star tree n=", n, " m=", m," strategy=", test_strategy, sep="")
      title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments. \n Star tree - setup ", setup, sep="")
    } else if (tree=="cat_binary"){
      name = paste(prefix, "_", "caterpillar", "_n=", n,  "_m=", m, sep="")
      subtitle = paste("Caterpillar tree n=", n, " m=", m," strategy=", test_strategy, sep="")
      title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments. \n Caterpillar tree", sep="")
    }
    
    # Plot & save
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
}

