library(foreach)
library(doParallel)
library(Rfast)  #rmvnorm
library(MASS)
setwd("src/star-tree")
source("tests.R")

#################
# Set variables #
#################

m = 20
n = 500
setup = 2
test_strategy = "half-and-half"  
# "half-and-half", "1-dependent" or "calculate-Y"

B = 3  # just for "1-dependent"
E = 1000
alphas = seq(0.01, 0.99, 0.01)

nr_exp = 500
save=TRUE



###############################################
# Compute empirical test sizes for all alphas #
###############################################

cores = 20#detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)


results <- foreach(nr = 1:nr_exp, .combine=rbind, .packages=c("MASS", "Rfast", "CombMSC", "RcppHelpers")) %dopar% {
  
                   
  if((nr%%10) == 0){
    print(nr)
  }
  
  # Generate n independent data sets depending on the setup
  if (setup==1){
    beta = rep(1,m)
    Sigma = beta %*% t(beta) + diag(rep(1,m))
  } 
  if (setup==2){
    #set.seed(nr+nr_exp)
    beta = c(10,10, rnorm((m-2),0,0.2))
    Sigma = beta %*% t(beta) + diag(rep(1/3,m))
  } 
  #set.seed(nr)
  X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma)
  
  
  # Call the test
  if (test_strategy=="half-and-half"){
    res = test_half_and_half(X, E=E, alphas=alphas)
    
  } else if (test_strategy=="1-dependent"){
    res = test_1_dependent(X, B=B, E=E, alphas=alphas)
    
  } else if (test_strategy=="calculate-Y"){
    res = test_calculate_Y(X, E=E, alphas=alphas)
    
  } else {
    print("ERROR")
  }
  as.numeric(res)
}


sizes = colMeans(results)
stopCluster(cl)



#########################
# Plot and save results #
#########################

subtitle = paste("Setup ", setup, ", (m,n) = (", m, ",", n, ")", sep="")
title = paste("Emprical test sizes vs. nominal test levels based on ", nr_exp, " experiments", sep="")

# Plot
if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "star-tree_setup", setup, "_", m, "_", n, ".pdf", sep="")
  name_rds = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "star-tree_setup", setup, "_", m, "_", n, ".rds", sep="")
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

