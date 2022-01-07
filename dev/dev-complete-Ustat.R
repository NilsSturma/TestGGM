library(MASS)
library(igraph)
library(TestGGM)

# General
n = 500
alpha = 0.05
E = 1000

# Tree
tree = "star_tree"  # Possible: "star_tree", "cat_binary"
m = 10
setup = 2

# Create tree
if (tree=="star_tree"){
  g = star_tree(m)
} else if (tree=="cat_binary"){
  g = cat_binary(m)
} 

# Plot tree
# pdf("dev-tree-plot.pdf") 
# plot(g)
# dev.off()

# Save all paths between all nodes in the tree
paths = get_paths(g)

# Collect the representations of the polynomial equality constraints that have to be tested
res = collect_indices(g, m)
ind_eq = res$ind_eq
p = dim(ind_eq)[1]
print(p)

# Generate n independent data sets
cov = cov_from_star_tree(g, paths, setup=setup, m=m)
X = mvrnorm(n, mu=rep(0,nrow(cov)), Sigma=cov)

test_complete_Ustat(X, ind_eq, E=1000)
