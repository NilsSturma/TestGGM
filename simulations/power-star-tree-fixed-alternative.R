library(foreach)
library(doParallel)
library(MASS) #mvrnorm
library(igraph)
library(TestGLTM)

source("simulations/utils.R") # TODO: add these functions to package

#################
# Set variables #
#################

# General
n_range = seq(50,1000, len=20)
E = 1000
nr_exp = 500
alpha = 0.05
m = 10  
setup = 2

beta_2 = c(rep(0,(m-2)),1,1)
h = 2


# Test strategy
test_strategy="symmetric" #"symmetric", "run-over", "U-stat"
B = 5  # just for test_strategy=="run-over" (5 works best for setup 1 after doing some experiments)
#beta = 0.001  # just for test_strategy=="two-step"
N = 2*n  # just for test_strategy=="U-stat"


  



# save
save=FALSE
name = paste(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_", "star-tree_setup=", setup, "_n=", n, "_m=", m, sep="")



###################################
# Create tree and collect indices #
###################################


g = star_tree(m)
plot(g)

res = collect_indices(g)
ind_eq = matrix(unlist(res[[1]]), ncol = 8, byrow = TRUE)
ind_ineq1 = matrix(unlist(res[[2]]), ncol = 6, byrow = TRUE)
ind_ineq2 = matrix(unlist(res[[3]]), ncol = 8, byrow = TRUE)

p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]



######################################
# Compute power for each alternative #
######################################

cores = 20  # detectCores()
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

results <- foreach(n = n_range, .combine=rbind, .packages=c("MASS", "TestGLTM", "igraph")) %dopar% {
  
  warnings()
  
  
  # Simulate power
  powers = rep(0, nr_exp)
  for (nr in 1:nr_exp){
    
    if((nr%%100) == 0){
      print(nr)
    }
    
    #Calculate covariance matric of alternative (depends on h)
    cov = cov_from_star_tree(g, setup=setup, m=m)
    cov = cov +  beta_2 %*% t(beta_2) * (h / sqrt(n))
    
    # Generate n indep datasets from the alternative
    X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
    
    # Call the test
    if (test_strategy=="run-over"){
      powers[nr] = test_m_dep(X, ind_eq, ind_ineq1, ind_ineq2, B=B, E=E, alphas=alpha)
    } else if (test_strategy=="symmetric"){
      powers[nr] = test_symmetric(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alpha)
    } else if (test_strategy=="U-stat"){
      powers[nr] = test_U_stat(X, ind_eq, ind_ineq1, ind_ineq2, N=N, E=E, alphas=alpha)
    }
  }
  simulated_power = mean(powers)
}

stopCluster(cl)


# Plot
if (save){
  # use "./img/name.png" to save in subdirectory
  name_pdf = paste("results/star-tree-general/", test_strategy, "/power-fixed-alternative/", name, ".pdf", sep="")
  name_rds = paste("results/star-tree-general/", test_strategy, "/power-fixed-alternative/", name, ".rds", sep="")
  saveRDS(results, file = name_rds) # read with readRDS()
  pdf(name_pdf) # create pdf file
}


title = paste("Emprical power for different n on a fixed alternative. \nStar tree - setup ", setup, sep="")
subtitle = paste("Local alternative = psi + b*t(b) + c*t(c) *", h, "/sqrt(n) with c=(rep(0,m-2),1,1).", sep="")
plot(n_range, results, 
     xlab="n", ylab="Emprical power", main=title, sub=subtitle,
     type="p", pch=1)
legend = c(paste("test-strategy = ", test_strategy, sep=""), 
           paste(nr_exp, " experiments", sep=""), 
           paste("m = ", m, sep=""))
legend("bottomright", legend =legend, bty = "n", lwd=0.1,
       cex = 1, lty = c(NA, NA, NA, NA))

if (save){
  dev.off() # close pdf file
}