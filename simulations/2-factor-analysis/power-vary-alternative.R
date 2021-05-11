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
n = 500
alpha = 0.05

# Test parameters
E = 1000
strategies = c("LR", "indep", "Ustat")  # Possible: "Ustat", "indep", "LR"
N = 5000

# Setup
m=20
setups= c("regular")  #do singular with H = seq(0.5,10,0.5) as well
nr_minors=10000

# Determine range of alternatives
H = seq(1,20,1)

# Parameter for simulations
nr_exp = 500
cores = 20
save = TRUE


for (strategy in strategies){
  for (setup in setups){
    
    # print some info
    print(paste("strategy = ",strategy ,sep=""))
    print(paste("setup = ", setup, sep=""))
    
    # initialize cores 
    cl <- makeCluster(cores, outfile = "")
    registerDoParallel(cl)
    
    # main computation
    results <- foreach(h = H, 
                       .combine=rbind, 
                       .errorhandling="stop",   
                       .packages=c("MASS", "TestGLTM", "stats")) %dopar% {
                         
       warnings()
       
       # Simulate power with nr_exp experiments. 
       powers = rep(0, nr_exp)
       for (nr in 1:nr_exp){
         
         # Print some info
         if((nr%%10) == 0){print(nr)}
         
         # Generate n indep datasets from the alternative
         cov = create_cov(setup, m, n, h=h)
         X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
         
         # Call test
         if (strategy=="LR"){
           res = factanal(X, 2)
         } else if (strategy=="indep"){
           res = test_indep_factors(X, nr_minors, E=E)
         } else if (strategy=="Ustat"){
           res = test_U_stat_factors(X, nr_minors, N=N, E=E)
         }
         powers[nr] = (res$PVAL <= alpha) # result: TRUE = rejected
       }
       simulated_power = mean(powers)
    }
    
    # Stop cluster
    stopCluster(cl)
    
    # Plot and save results
    name = paste(setup, "_n=", n, "_m=", m, sep="")
    subtitle = paste("Based on ", nr_exp, " experiments.", sep="")
    if (save){
      name_pdf = paste("./results/2-factor/power/", strategy, "/vary-alternative_", name, ".pdf", sep="")
      name_rds = paste("./results/2-factor/power/", strategy, "/vary-alternative_", name, ".rds", sep="")
      saveRDS(sizes, file = name_rds)
      pdf(name_pdf)
    }
    plot(H, results, 
         xlab="h", ylab="power", sub=subtitle,
         type="p", pch=1, ylim=c(0,1))
    abline(coef = c(0,1))
    if (save){dev.off()}
  }
}


