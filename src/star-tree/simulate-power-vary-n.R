library(foreach)
library(doParallel)
library(MASS)
source("tests.R")

#################
# Set variables #
#################

m = 10
setup = 1

B = 3 
E = 1000
alpha = 0.05
n_range = seq(100,1000,50)
nr_exp = 500

save=TRUE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "fixed-alternative_alpha=", alpha, "_vary-n", sep="")


###################
# Fix alternative #
###################

# Generate covariance matrix depending on the setup
if (setup==1){
  beta = rep(1,m)
  Sigma = beta %*% t(beta) + diag(rep(1,m))
} 
if (setup==2){
  beta = c(10,10, rnorm((m-2),0,0.2))
  Sigma = beta %*% t(beta) + diag(rep(1/3,m))
} 

# Create noise
a = -0.5
b = 0.5
shift_unif = function(x,a,b){a + (b-a) * x}

A = matrix(shift_unif(runif(m**2), a, b), ncol=m)   # elements in (a,b)
noise = t(A) %*% A # alwas pos. semidefinit

# add noise to final covariance matric
Sigma = Sigma + noise


#############################
# Simulate power for each n #
#############################

cores = detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)


results <- foreach(n = n_range, .combine=rbind, .packages=c("MASS", "CombMSC")) %dopar% {
  
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

title = paste("Emprical power for different n. Each based on ", nr_exp, " experiments. Alpha = ", alpha, sep="")
subtitle = paste("Fixed alternative: Setup ", setup, " with noise Unif[", a, ",", b, "].", sep="")

if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste(name, ".pdf", sep="")
  name_rds = paste(name, ".rds", sep="")
  saveRDS(sizes, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}

plot(n_range, results, 
     xlab="n", ylab="Emprical power", main=title, sub=subtitle,
     type="p", pch=1)
abline(coef = c(0,1))

if (save){
  dev.off() # close pdf file
}