library(foreach)
library(doParallel)
library(MASS)
source("tests.R")

#################
# Set variables #
#################

m = 10
n = 500 # fixed
setup = 1

B = 3 
E = 1000
alpha = 0.05 # fixed
nr_exp = 500

# alternative
beta_2 = c(1,1, rep(0,(m-2)))
H = seq(1, 20, 2)

save=TRUE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_n=", n, "_alpha=", alpha, "_vary-alternative", sep="")


#############################
# Simulate power for each n #
#############################

cores = 20#detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)


results <- foreach(h = H, .combine=rbind, .packages=c("MASS", "CombMSC")) %dopar% {
  
  # Calculate covariance matrix for alternative
  if (setup==1){
    beta = rep(1,m)
    Sigma = beta %*% t(beta) + diag(rep(1,m))
  } 
  if (setup==2){
    beta = c(10,10, rnorm((m-2),0,0.2))
    Sigma = beta %*% t(beta) + diag(rep(1/3,m))
  } 
  
  
  # add noise to final covariance matrix
  Sigma = Sigma + beta_2 %*% t(beta_2) * h / sqrt(n) 
  
  # Simulate power
  powers = rep(0, nr_exp)
  for (nr in 1:nr_exp){
    
    if((nr%%10) == 0){
      print(nr)
    }
    
    # Generate n independent data sets from the alternative
    X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma)
    
    # Call the test
    res = test_1_dependent(X, B=B, E=E, alphas=alpha)
    powers[nr] = res
  }
  simulated_power = mean(powers)
}

stopCluster(cl)



#########################
# Plot and save results #
#########################

title = paste("Emprical power for different deviations of the star tree.\nEach based on ", nr_exp, " experiments.", sep="")
subtitle = paste("Deviations from setup ", setup, ", n = ", n, ", alpha = ", alpha, ".", sep="")

if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste(name, ".pdf", sep="")
  name_rds = paste(name, ".rds", sep="")
  saveRDS(results, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}

plot(H, results, 
     xlab="deviation", ylab="Emprical power", main=title, sub=subtitle,
     type="p", pch=1)

if (save){
  dev.off() # close pdf file
}
