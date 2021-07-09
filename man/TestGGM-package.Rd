\name{TestGGM-package}
\alias{TestGGM-package}
\alias{TestGGM}
\docType{package}
\title{
  Testing Gaussian Graphical Models
}
\description{
  Tests the goodness-of-fit of Gaussian latent tree models and the two-factor analysis model. Three algebraic testing methodologies are implemented as well as the likelihood ratio test using the EM-algorithm. Supports all Gaussian latent tree models where the observed variables correspond to the leaves of the tree as well as the factor analysis model.
}
\author{
Nils Sturma <nils.sturma@tum.de>
}
%\details{}
%\keyword{}
\references{
TO BE WRITTEN
}
\examples{
  \dontrun{
      # Create tree
      vertices <- data.frame(name=seq(1,8), type=c(rep(1,5), rep(2,3))) # 1=observed, 2=latent
      edges <- data.frame(from=c(1,2,3,4,5,6,7), to=c(8,8,6,6,7,7,8))
      tree <- graph_from_data_frame(edges, directed=FALSE, vertices=vertices)
      
      # Sample data from tree
      V(tree)$var = rep(1,8)
      E(tree)$corr = rep(0.7,7)
      X = sample_from_tree(tree, m=5, n=500)
      
      # Goodness of fit test
      gltmtest(X, tree, m=5, test_strategy="grouping")
  }
}
