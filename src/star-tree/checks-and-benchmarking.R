library(MASS)
library(CombMSC)
library(microbenchmark)
library(RcppHelpers)




################
# benchmarking #
################

# generate data
n = 500
m = 10
E = 1000
B = 3

beta = rep(1,m)
Sigma = beta %*% t(beta) + diag(rep(1,m))
X = mvrnorm(n, mu=rep(0,m), Sigma=Sigma)

omega = floor((n-1)/B)
X1 = X[1:(n-1),]
X2 = X[2:n,]

# Create indices
sub_sets = subsets(m,4,1:m)
nr_cols = 2*nrow(sub_sets)
indices = matrix(0, nrow=nr_cols, 4)
indices[1:(nr_cols/2),] = sub_sets
indices[((nr_cols/2)+1):(2*(nr_cols/2)),] = sub_sets[,c(1,3,2,4)]
mode(indices) = "integer"

# Compute Y
Y = calculate_Y(indices, X1, X2) # each column is one polynom
Y_mean = colmeans(Y)
Y_centered = transpose(transpose(Y) - Y_mean) # Centering: Y_i = (Y_i - Y_mean)

# Compute the diagonal of the batched mean estimator cov(Y)
cov_Y_diag = rep(0, length(Y_mean))
for (b in 1:omega){
  L = seq(1+(b-1)*B, b*B)
  cov_Y_diag = cov_Y_diag + colsums(Y_centered[L,])**2
}
cov_Y_diag = cov_Y_diag / (B*omega)

# Vector for standardizing
standardizer = cov_Y_diag**(-1/2)

# Test statistic
test_stat = sqrt(n-1) * max(abs(standardizer * Y_mean))


# Bootstrapping 
bootstrapR = function(E, B, omega, standardizer, Y_centered){
  results = rep(0, E)
  for (i in 1:E){
    epsilons = Rnorm(omega, m=0, s=1)
    sum = rep(0, length(Y_mean))
    for (b in 1:omega){
      L = seq(1+(b-1)*B, b*B)
      sum = sum + epsilons[b] * colsums(Y_centered[L,])
    }
    results[i] = (1/sqrt(B*omega)) * max(abs(standardizer * sum))
  }
  return(results)
}

resR = bootstrapR(E, B, omega, standardizer, Y_centered)
res = bootstrap_dependent(E, B, omega, standardizer, Y_centered)
  


microbenchmark(resR = bootstrapR(E, B, omega, standardizer, Y_centered),
               res = bootstrap_dependent(E, B, omega, standardizer, Y_centered),
               times=10L)
