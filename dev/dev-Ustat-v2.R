library(MASS)
library(igraph)
library(TestGGM)

# General
n = 500
alpha = 0.05
E = 1000

# Tree
tree = "star_tree"  # Possible: "star_tree", "cat_binary"
m = 20
setup = 1

# Create tree
if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="cat_binary"){
  g = cat_binary(m)
} 

# Save all paths between all nodes in the tree
paths = get_paths(g)

# Collect the representations of the polynomial equality constraints that have to be tested
res = collect_indices(g, m)
ind_eq = res$ind_eq
ind_ineq1 = res$ind_ineq1
ind_ineq2 = res$ind_ineq2
p = dim(ind_eq)[1] + dim(ind_ineq1)[1] + dim(ind_ineq2)[1]
print(p)

# Determine computational budget
N = floor(log(p) * n^(2/3))
print(N)

# Covariance matrix
cov = cov_from_star_tree(g, paths, setup=setup, m=m)

# Apply test
X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)
test_U_stat_v2(X, ind_eq, ind_ineq1, ind_ineq2, N=125, E=E)
