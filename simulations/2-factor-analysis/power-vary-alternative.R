library(MASS)
library(stats)
library(rlist)
library(foreach)
library(TestGLTM)
library(doParallel)
setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")
source("./simulations/2-factor-analysis/utils.R")






################
## Parameters ##
################


# Sample size
n = 500
alpha = 0.05

# Test parameters
E = 1000
strategies = c("indep")  # Possible: "Ustat", "indep", "LR"
N = 5000

# Setup
m=20
setups= c("regular", "singular")  #do singular with H = seq(0.5,10,0.5) as well
nr_minors=10000
randomized=TRUE

# Determine range of alternatives
H = seq(1,60,1)

# Parameter for simulations
nr_exp = 500
cores = 20
save = TRUE


create_minors <- function(m, randomized=FALSE, nr_minors=10000){
  
  if (randomized){
    # choose sets of indices
    A = matrix(unlist(random_combs(m,6,nr_minors)[[1]]), 
               ncol = 6, byrow = TRUE)
    # shuffle each row
    A = t(apply(A, 1, sample))
    ind_minors = cbind(t(apply(A[,1:3],1,sort)), t(apply(A[,4:6],1,sort)))
    
  } else {
    sub_sets = subsets(m,6,1:m)
    L = list()
    for (s in 1:dim(sub_sets)[1]){
      set = sub_sets[s,]
      for (i in 2:5){
        for (j in (i+1):6){
          C = set[c(1,i,j)]
          B = setdiff(set, C)
          L = list.append(L, c(C,B))
        }
      }
    }
    ind_minors = matrix(unlist(L), ncol = 6, byrow = TRUE)
  }
  
  
  return(ind_minors)
}
ind_minors = create_minors(m, randomized=randomized, nr_minors=nr_minors)




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
         cov = create_cov(setup, m, n, h)
         X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
         
         # Call test
         if (strategy=="LR"){
           res = factanal(X, 2)
         } else if (strategy=="indep"){
           res = test_indep_factors(X, ind_minors, E=E)
         } else if (strategy=="Ustat"){
           res = test_U_stat_factors(X, ind_minors, N=N, E=E)
         }
         powers[nr] = (res$PVAL <= alpha) # result: TRUE = rejected
       }
       simulated_power = mean(powers)
    }
    
    # Stop cluster
    stopCluster(cl)
    
    # Plot and save results
    name = paste("hlarge_", setup, "_n=", n, "_m=", m, sep="")
    subtitle = paste("Based on ", nr_exp, " experiments.", sep="")
    if (save){
      name_pdf = paste("./results/2-factor/power/", strategy, "/vary-alternative_", name, ".pdf", sep="")
      name_rds = paste("./results/2-factor/power/", strategy, "/vary-alternative_", name, ".rds", sep="")
      saveRDS(results, file = name_rds)
      pdf(name_pdf)
    }
    plot(H, results, 
         xlab="h", ylab="Empirical power", sub=subtitle,
         type="p", pch=1, ylim=c(0,1))
    if (save){dev.off()}
  }
}


