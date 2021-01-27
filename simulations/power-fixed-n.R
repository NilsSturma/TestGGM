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
n = 500
E = 1000
nr_exp = 500
alpha = 0.05

# Test strategy
strategies = c("grouping")  # Possible: "grouping", "run-over", "U-stat", "LR", "grouping-cov", "run-over-cov"
B = 5  # only relevant if test_strategy=="run-over" 
N = 5000  # only relevant if test_strategy=="U-stat"

# Tree
tree = "star_tree"  # Possible: "star_tree", "cat_binary"
m = 20  
setup = 2  # only relevant if tree=="star_tree"

# Determine the alternatives
beta_2 = c(rep(0,(m-2)),1,1)
H = seq(0.5,10,len=20)

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



######################################
# Compute power for each alternative #
######################################

for (test_strategy in strategies){
  
  print(test_strategy)
  
  # Initialize the cluster
  cores = 20  # detectCores()
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  # Main loop to compute the empirical power. For each h this is done one a separate compute node of the cluster.
  results <- foreach(h = H, 
                     .combine=rbind, 
                     .errorhandling="remove",
                     .packages=c("MASS", "TestGLTM", "igraph", "stats")) %dopar% {
    
    warnings()
    
    # Simulate power with nr_exp experiments. 
    powers = rep(0, nr_exp)
    for (nr in 1:nr_exp){
      
      # Print some info
      if((nr%%20) == 0){print(nr)}
      
      # Calculate covariance matric of alternative (depends on h)
      if (tree=="star_tree"){
        cov = cov_from_star_tree(g, paths, setup=setup, m=m)
      } else if (tree=="cat_binary"){
        V(g)$var = rep(1,38)
        E(g)$corr = rep(0.7,37)
        cov = cov_from_graph(g, m, paths)
      }
      cov = cov +  beta_2 %*% t(beta_2) * (h / sqrt(n))
      
      # Generate n indep datasets from the alternative
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
      } else if (test_strategy=="grouping-cov"){
        res = test_grouping_compute_cov(X, ind_eq, E=E) # only possible for equalities
      } else if (test_strategy=="run-over-cov"){
        res = test_run_over_compute_cov(X, ind_eq, E=E) # only possible for equalities
      }
      
      # Rejected?
      powers[nr] = res$PVAL <= alpha # result: TRUE = rejected
    }
    
    simulated_power = mean(powers)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Plot and save result
  # Create meta information for the plot and the files to be saved
  if (tree=="star_tree"){
    name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "setup=", setup, "_n=", n, "_m=", m, sep="")
    title = paste("Emprical power for different local alternatives based on ",  nr_exp, 
                  " experiments. \n Star tree - setup ", setup, ", strategy=", test_strategy, sep="")
    subtitle = "Local alternative = psi + b*t(b) + c*t(c) *h/sqrt(n) with c=c(rep(0,(m-2)),1,1)."
  } else if (tree=="cat_binary"){
    name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "caterpillar", "_n=", n,  sep="")
    title = paste("Emprical power for different local alternatives based on ",  nr_exp, 
                  " experiments. \n Caterpillar tree, strategy=", test_strategy, sep="")
    subtitle = "Local alternative = Sigma + c*t(c) *h/sqrt(n) with c=c(rep(0,(m-2)),1,1)."
  }
  
  # Plot
  if (save){
    # use "./img/name.png" to save in subdirectory
    name_pdf = paste("./results/", tree, "/", test_strategy, "/power-fixed-n/", name, ".pdf", sep="")
    name_rds = paste("./results/", tree, "/", test_strategy, "/power-fixed-n/", name, ".rds", sep="")
    saveRDS(results, file = name_rds) # read with readRDS()
    pdf(name_pdf) # create pdf file
  }
  plot(H, results, 
       xlab="h", ylab="Emprical power", main=title, sub=subtitle,
       type="p", pch=1)
  if (save){
    dev.off() # close pdf file
  }
}