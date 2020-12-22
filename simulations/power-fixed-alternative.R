library(foreach)
library(doParallel)
library(MASS) #mvrnorm
library(igraph)
library(TestGLTM)

setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")
source("simulations/utils.R") # TODO: add these functions to package

#################
# Set variables #
#################

# General
n_range = seq(250,1200, 50)
E = 1000
nr_exp = 500
alpha = 0.05

# Tree
tree = "cat_binary"  # "star_tree", "cat_binary"
m = 20  
setup = 2  # (star_tree)


beta_2 = c(rep(0,(m-2)),1,1)
h = 15


# Test strategy
test_strategy="run-over"  # "grouping", "run-over", "U-stat", "LR"
B = 5  # just for test_strategy=="run-over" (5 works best for setup 1 after doing some experiments)
N = 5000  # just for test_strategy=="U-stat"



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



######################################
# Compute power for each alternative #
######################################




results = rep(0, length(n_range))
for (i in (1:length(n_range))){
  
  cores = 20  # detectCores()
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  n = n_range[i]
  print(paste("n=",n ,sep=""))
  
  # Simulate power
  powers <- foreach(nr = 1:nr_exp, 
                    .combine=rbind, 
                    .errorhandling="remove",
                    .packages=c("MASS", "TestGLTM", "igraph", "stats")) %dopar% {
    
    #warnings()
    if((nr%%20) == 0){
      print(nr)
    }
    
    
    # Calculate covariance matric of alternative (depends on h)
    if (tree=="star_tree"){
      cov = cov_from_star_tree(g, setup=setup, m=m)
    } else if (tree=="cat_binary"){
      V(g)$var = rep(1,38)
      E(g)$corr = rep(0.7,37)
      cov = cov_from_graph(g)
    }
    
    cov = cov +  beta_2 %*% t(beta_2) * (h / sqrt(n))
    
    # Generate n indep datasets from the alternative
    X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
    
    # Call the test
    if (test_strategy=="LR"){
      if (tree=="star_tree"){
        res = factanal(X, 1)
        result = res[["PVAL"]] <= alpha # result: TRUE = rejected
      } else if (tree=="cat_binary"){
        result = LR_test(X,g) <= alpha # result: TRUE = rejected
      }
    } else if (test_strategy=="grouping"){
      result = test_grouping(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alpha)
    } else if (test_strategy=="run-over"){
      result = test_run_over(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E, alphas=alpha)
    } else if (test_strategy=="U-stat"){
      result = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E, alphas=alpha)
    } else if (test_strategy=="U-stat-deg"){
      result = test_U_stat_degenerate(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E, alphas=alpha)
    }
    result = as.numeric(result)
  }
  results[i] = mean(powers)
  print(results[i])
  stopCluster(cl)
}




#########################
# Plot and save results #
#########################


if (tree=="star_tree"){
  name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "setup=", setup, "_m=", m, sep="")
  title = paste("Emprical power for different n with fixed alternative based on ",  nr_exp, 
                " experiments. \n Star tree - setup ", setup, ", strategy=", test_strategy, sep="")
  subtitle = paste("Local alternative = psi + b*t(b) + c*t(c) *", h ,"/sqrt(n) with c=c(rep(0,(m-2)),1,1).", sep="")
} else if (tree=="cat_binary"){
  name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "caterpillar", sep="")
  title = paste("Emprical power for different n with fixed alternative based on ",  nr_exp, 
                " experiments. \n Caterpillar tree, strategy=", test_strategy, sep="")
  subtitle = paste("Local alternative = Sigma + c*t(c) *15/sqrt(n) with c=c(rep(0,(m-2)),1,1).", sep="")
}



# Plot
if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste("./results/", tree, "/", test_strategy, "/power-fixed-alternative/", name, ".pdf", sep="")
  name_rds = paste("./results/", tree, "/", test_strategy, "/power-fixed-alternative/", name, ".rds", sep="")
  saveRDS(results, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}




plot(n_range, results, 
     xlab="h", ylab="Emprical power", main=title, sub=subtitle,
     type="p", pch=1)
#legend = c(paste("test-strategy = ", test_strategy, sep=""), 
#           paste(nr_exp, " experiments", sep=""), 
#           paste("n = ", n, sep=""), 
#           paste("m = ", m, sep=""))
#legend("bottomright", legend =legend, bty = "n", lwd=0.1,
#       cex = 1, lty = c(NA, NA, NA, NA))

if (save){
  dev.off() # close pdf file
}

