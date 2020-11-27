library(TestGLTM)
library(stats)
library(MASS)

n = 250
nr_exp = 500
alphas = seq(0.01, 0.99, 0.01)

m = 20  
setup = 2
g = star_tree(m)

cores = 20  # detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

results <- foreach(nr = 1:nr_exp, .combine=rbind, .packages=c("MASS", "TestGLTM", "igraph")) %dopar% {
  
  if((nr%%10) == 0){
    print(nr)
  }
  warnings()
  
  # Generate n independent data sets depending on setup
  cov = cov_from_star_tree(g, setup=setup, m=m)
  X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
  
  # Call the test
  res = factanal(X, 1)
  is_rejected = res[["PVAL"]] <= alphas
  is_rejected = as.numeric(is_rejected)
}

sizes = colMeans(results)
stopCluster(cl)

plot(alphas, sizes, 
     xlab="Nominal level", ylab="Emprical test size", #main=title, sub=subtitle,
     type="p", pch=1)
abline(coef = c(0,1))