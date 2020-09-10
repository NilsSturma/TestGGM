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
n_range = seq(100,1000,50)
nr_exp = 500

deviations = seq(0.1, 1.5, 0.1)

save=TRUE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_n=", n, "_alpha=", alpha, "_vary-alternative", sep="")

shift_unif = function(x,a,b){a + (b-a) * x}


#############################
# Simulate power for each n #
#############################

cores = detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)


results <- foreach(d = deviations, .combine=rbind, .packages=c("MASS", "CombMSC")) %dopar% {
  
  # Calculate covariance matrix for alternative
  if (setup==1){
    beta = rep(1,m)
    Sigma = beta %*% t(beta) + diag(rep(1,m))
  } 
  if (setup==2){
    beta = c(10,10, rnorm((m-2),0,0.2))
    Sigma = beta %*% t(beta) + diag(rep(1/3,m))
  } 
  
  # Create noise
  a = -d
  b = d
  A = matrix(shift_unif(runif(m**2), a, b), ncol=m)   # elements in (a,b)
  noise = t(A) %*% A # alwas pos. semidefinit
  
  # add noise to final covariance matrix
  Sigma = Sigma + noise
  
  # Simulate power
  powers = rep(0, nr_exp)
  for (nr in 1:nr_exp){
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

title = paste("Emprical power for different deviations of the star tree. Each based on ", nr_exp, " experiments.", sep="")
subtitle = paste("Deviations from setup ", setup, ", n = ", n, ", alpha = ", alpha, ".", sep="")

if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste(name, ".pdf", sep="")
  name_rds = paste(name, ".rds", sep="")
  saveRDS(results, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}

plot(n_range, results, 
     xlab="deviation", ylab="Emprical power", main=title, sub=subtitle,
     type="p", pch=1)
abline(coef = c(0,1))

if (save){
  dev.off() # close pdf file
}
