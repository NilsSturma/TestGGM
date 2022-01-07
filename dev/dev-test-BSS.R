library(MASS)
library(igraph)
library(TestGGM)

# General
n = 5000
alphas = alphas = seq(0.01, 0.99, 0.01)
betas =  alphas/10
E = 1000

# Tree
tree = "star_tree"  # Possible: "star_tree", "cat_binary"
m = 8
setup = 2

# Create tree
if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="cat_binary"){
  g = cat_binary(m)
} 

# Plot tree
pdf("dev-tree-plot.pdf") 
plot(g)
dev.off()

# Save all paths between all nodes in the tree
paths = get_paths(g)

# Collect the representations of the polynomials that have to be tested
res = collect_indices(g, m)
ind_eq = res$ind_eq
ind_ineq1 = res$ind_ineq1
ind_ineq2 = res$ind_ineq2
p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
print(p)

# Generate n independent data sets depending on setup
if (tree=="star_tree"){
    cov = cov_from_star_tree(g, paths, setup=setup, m=m)
} else if (tree=="cat_binary"){
    V(g)$var = rep(1,(m+(m-2))) #parameters & starting values of EM algo
    E(g)$corr = rep(0.7,(m+(m-3)))
    cov = cov_from_graph(g, m, paths)
}
X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)

begin = Sys.time()
test_grouping_BSS(X, ind_eq, ind_ineq1, ind_ineq2, E=E, alphas=alphas, betas=betas)
end = Sys.time()
print(end-begin)

#########
## BSS ##
#########
N = findn(nrow(X),4)
indices = matrix(1:N, ncol=4)
H = calculate_H(X, indices, ind_eq, ind_ineq1, ind_ineq2)

n = dim(H)[1]
p = dim(H)[2]  # total nr of constraints
p_eq = dim(ind_eq)[1]  # nr of equality constraints

# Mean and centering
H_mean = Rfast::colmeans(H)
H_centered = Rfast::transpose(Rfast::transpose(H) - H_mean) # Centering: H_i = (H_i - H_mean)

# Diagonal of the sample covariance of H
cov_H_diag = Rfast::colsums(H_centered**2) / n

# Vector for standardizing
standardizer = cov_H_diag**(-1/2)

# Test statistic
marginal_stats = sqrt(n) * H_mean
marginal_stats[1:p_eq] = abs(marginal_stats[1:p_eq])
test_stat =  max(standardizer * marginal_stats)

# Bootstrapping - first step  (only inequalities)
W = bootstrap(E, H_centered[,(p_eq+1):p])

# Calculate c_beta
W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer[(p_eq+1):p])
bootstrap_res = Rfast::rowMaxs(W_standardized, value = TRUE)
c_beta = as.numeric(quantile(bootstrap_res, probs=1-beta))


# Calculate nuisance parameter lambda (zero for all equalities)
lambda = H_mean[(p_eq+1):p] + (cov_H_diag[(p_eq+1):p]**(1/2)) * c_beta/sqrt(n)
lambda = sapply(lambda, min_zero)

 # Bootstrapping - second step
W = bootstrap(E, H_centered)
W[,1:p_eq] = abs(W[,1:p_eq])
W[,(p_eq+1):p] = Rfast::transpose(Rfast::transpose(W[,(p_eq+1):p]) + lambda * sqrt(n))
W_standardized = Rfast::transpose(Rfast::transpose(W) * standardizer)
maxima = Rfast::rowMaxs(W_standardized, value = TRUE)
critical_value = as.numeric(quantile(maxima, probs=1-alpha+beta))

# Reject?
test_stat > critical_value