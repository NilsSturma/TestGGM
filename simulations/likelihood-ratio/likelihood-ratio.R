library(TestGLTM)
library(MASS)  # mvrnorm
library(matrixcalc)  # is.positive.definite
source("simulations/utils.R")

n = 250
m = 10  
setup = 1
g = star_tree(m)
plot(g)
cov = cov_from_star_tree(g, setup=setup, m=m)
X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)


########################
# function definitions #
########################

vec2mat <-function(v){
  m = (-1 + sqrt(1+8*55)) / 2
  mat1 = matrix(, m, m)
  mat1[lower.tri(mat1, diag=TRUE)] = v
  mat2 = t(mat1)
  mat2[lower.tri(mat2, diag=TRUE)] = v
  return(mat2)
}


mat2vec = function(M){
  return(as.vector(M[lower.tri(M, diag=TRUE)]))
}


sample_cov <- function(X){
  n = dim(X)[1]
  (1/n) * t(X) %*% X
}


loglik <- function(X, S){
  if (!is.positive.definite((S))){
    return(-1000000) # very low value
  }
  n = dim(X)[1]
  res = - (n/2) * log(det(S)) - (1/2) * sum(diag( t(X)%*%X %*% solve(S)))
  # Skipped constants in likelihood function
  return(res)
}

eval_loglik <- function(v){
  S = vec2mat(v)
  return(-loglik(X,S))  # to be minimized
}


eq_constr <- function(S, ind_eq){
  constr = vector("numeric",length = nrow(ind_eq))
  for (j in 1:nrow(ind_eq)) {
    constr[j] = S[ind_eq[j,1], ind_eq[j,2]] * S[ind_eq[j,3],ind_eq[j,4]] - S[ind_eq[j,5], ind_eq[j,6]] * S[ind_eq[j,7],ind_eq[j,8]]
  }
  return(constr)
}


eval_eq_constr <- function(v){
  S = vec2mat(v)
  return(eq_constr(S,ind_eq))
}



# The sample covariance maximizes the loglik
S = sample_cov(X)
global_max = loglik(X,S)



res = collect_indices(g)
ind_eq = matrix(unlist(res[[1]]), ncol = 8, byrow = TRUE)
ind_ineq1 = matrix(unlist(res[[2]]), ncol = 6, byrow = TRUE)
ind_ineq2 = matrix(unlist(res[[3]]), ncol = 8, byrow = TRUE)


library(nloptr)
# start value
x0 = mat2vec(sample_cov(X))
lower = mat2vec(matrix(-3, m, m))
upper = mat2vec(matrix(3, m, m))

#res = isres(x0, fn=eval_loglik, lower=lower, upper=upper, heq=eval_eq_constr)


opt <- nloptr(x0=x0, eval_f=eval_loglik, lb = lower, ub = upper,  eval_g_eq = eval_eq_constr, 
              opts = list(algorithm="NLOPT_GN_ISRES",maxeval=100000))




######################
# Example 
######################



F <- function(x,y,A){  # made-up function
  # minimize scaled distance between points x and y
  sum((A[1]*x-A[2]*y)^2)
}
Gc <- function(x,y,A) return(sum(y/3) - sum(x*y))
Hc <- function(x,y,A) return(1-sum(x))

library(nloptr)
y= c(0,1,0)
A= c(10,1)
opt <- nloptr(x0=rep(1/3,3), eval_f=F, lb = rep(0.05,3), ub = rep(1,3), 
              eval_g_ineq = Gc, eval_g_eq = Hc, 
              opts = list(algorithm="NLOPT_GN_ISRES",maxeval=100000), y=y, A=A)






