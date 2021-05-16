# Covariance in different setups
create_cov <- function(setup="regular", m, n, h=0){
  if (setup=="regular"){
    beta = matrix(stats::rnorm(2*m),m,2)
    Sigma = rep(1,m)
    cov = diag(Sigma) + beta %*% t(beta)
  } else {
    beta_1 = rep(1,m)#stats::rnorm(m,0,1)
    beta_2 = c(10,10, stats::rnorm((m-2),0,0.2))
    Sigma = rep(1/3,m)
    cov = diag(Sigma) + beta_1 %*% t(beta_1) + beta_2 %*% t(beta_2)
  }
  if (h!=0){ #create covariance for local alternative
    beta_3 = c(rep(0,(m-2)),1,1)
    cov = cov +  beta_3 %*% t(beta_3) * (h / sqrt(n))
  }
  return(cov)
}