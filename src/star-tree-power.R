library(foreach)
library(doParallel)
library(MASS)
source("star-tree-tests.R")

# fix alternative and alpha, vary n

#################
# Set variables #
#################

m = 10
setup = 1

B = 3
E = 1000
alpha = 0.05

nr_exp = 500
save=FALSE

seq(100, 1000, 50)

#################################################
# Generate covariance matrix sepending on setup #
#################################################
if (setup==1){
  beta = rep(1,m)
  Sigma = beta %*% t(beta) + diag(rep(1,m))
} 
if (setup==2){
  beta = c(10,10, rnorm((m-2),0,0.2))
  Sigma = beta %*% t(beta) + diag(rep(1/3,m))
} 


# Create noise
shift_unif = function(x,a,b){a + (b-a) * x}

n = 4
a = -0.5
b = 0.5
A = matrix(shift_unif(runif(n**2), a, b), ncol=n)   # elements in (-1,1)
noise = t(A) %*% A # alwas pos. semidefinit

X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma)