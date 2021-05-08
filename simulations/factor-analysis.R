library(MASS)
library(stats)
library(foreach)
library(TestGLTM)
library(doParallel)
setwd("/dss/dsshome1/lxc0D/ge73wex3/master-thesis-tests")


################
## Parameters ##
################


# Sample size
n = 500
alphas = seq(0.01, 0.99, 0.01)

# Test parameters
E = 1000
strategy = "Ustat"  # Possible: "Ustat", "indep", "LR"
N = 5000

# Setup
m=20
setup="singular"
nr_minors=10000

# Parameter for simulations
nr_exp = 500
cores = 20
save = FALSE



####################################
## Covariance in different setups ##
####################################
create_cov <- function(setup="regular", m=20){
  if (setup=="regular"){
    factors = 2
    beta = matrix(stats::rnorm(2*m),m,2)
    Sigma = rep(1,m)
    cov = diag(Sigma) + beta %*% t(beta)
  } else {
    beta_1 = rep(1,m)#stats::rnorm(m,0,1)
    beta_2 = c(10,10, stats::rnorm((m-2),0,0.2))
    Sigma = rep(1/3,m)
    cov = diag(Sigma) + beta_1 %*% t(beta_1) + beta_2 %*% t(beta_2)
  }
  return(cov)
}




#################
## Simulations ##
#################
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

results <- foreach(nr = 1:nr_exp, 
                   .combine=rbind, 
                   .errorhandling="stop",   
                   .packages=c("MASS", "TestGLTM", "stats")) %dopar% {
  
  # Print some info
  if((nr%%20) == 0){print(nr)}
  warnings()
  
  # sample parameters
  cov = create_cov(setup, m)
  
  # sample data
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  
  if (strategy=="LR"){
    res = factanal(X, 2)
  } else if (strategy=="indep"){
    res = test_indep_factors(X, nr_minors, E=E)
  } else if (strategy=="Ustat"){
    res = test_U_stat_factors(X, nr_minors, N=N, E=E)
  }
  res = (res$PVAL <= alphas) # result: TRUE = rejected
}

# Stop cluster
stopCluster(cl)

# Compute empirical sizes
sizes = colMeans(results)


# plot
plot(alphas, sizes, 
     xlab="Nominal level", ylab="Emprical test size",
     type="p", pch=1, ylim=c(0,1))
abline(coef = c(0,1))
