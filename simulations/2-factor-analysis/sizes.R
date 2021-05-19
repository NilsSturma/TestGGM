library(MASS)
library(stats)
library(foreach)
library(TestGLTM)
library(doParallel)
source("utils.R")
setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")


################
## Parameters ##
################


# Sample size
n_range = c(500)
alphas = seq(0.01, 0.99, 0.01)

# Test parameters
E = 1000
strategies = c("LR")  # Possible: "Ustat", "indep", "LR"
N = 5000

# Setup
m=200
setups= c("regular")  # Possible: regular, singular
nr_minors=10000

# Parameter for simulations
nr_exp = 500
cores = 20
save = FALSE






#################
## Simulations ##
#################

for (strategy in strategies){
  for (setup in setups){
    for (n in n_range){
      
      # print some info
      print(paste("strategy = ",strategy ,sep=""))
      print(paste("setup = ", setup, sep=""))
      print(paste("n = ",n ,sep=""))
      
      # initialize cores 
      cl <- makeCluster(cores, outfile = "")
      registerDoParallel(cl)
      
      # main computation
      results <- foreach(nr = 1:nr_exp, 
                         .combine=rbind, 
                         .errorhandling="stop",   
                         .packages=c("MASS", "TestGLTM", "stats")) %dopar% {
         # Print some info
         if((nr%%20) == 0){print(nr)}
         warnings()
         
         # Sample parameters and data 
         cov = create_cov(setup, m, n, h=0)
         X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
         
         # Call test
         if (strategy=="LR"){
           res = factanal(X, 2)
         } else if (strategy=="indep"){
           res = test_indep_factors(X, nr_minors, E=E)
         } else if (strategy=="Ustat"){
           res = test_U_stat_factors(X, nr_minors, N=N, E=E)
         }
         res = (res$PVAL <= alphas) # TRUE = rejected
      }
      
      # Stop cluster
      stopCluster(cl)
      
      # Compute empirical sizes
      sizes = colMeans(results)
      
      # Plot and save results
      name = paste(setup, "_n=", n, "_m=", m, sep="")
      subtitle = paste("Based on ", nr_exp, " experiments.", sep="")
      if (save){
        name_pdf = paste("./results/2-factor/size/", strategy, "/", name, ".pdf", sep="")
        name_rds = paste("./results/2-factor/size/", strategy, "/", name, ".rds", sep="")
        saveRDS(sizes, file = name_rds)
        pdf(name_pdf)
      }
      plot(alphas, sizes, 
           xlab="Nominal level", ylab="Empirical test size", sub=subtitle,
           type="p", pch=1, ylim=c(0,1))
      abline(coef = c(0,1))
      if (save){dev.off()}
    }
  }
}