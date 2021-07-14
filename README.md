## TestGGM
R package for testing the goodness-of-fit of Gaussian graphical models. Supports two types of models:
- Gaussian latent tree models
- Two-factor analysis model

Install via `devtools::install_github("NilsSturma/TestGGM")`. The package implements the tests described in TO BE WRITTEN.

Most important files:
- `R/tests.R`: algebraic tests for (arbitrary) Gaussian latent tree models 
- `R/LR_EM.R`: likelihood-ratio test for (arbitrary) Gaussian latent tree models using the EM algorithm
- `R/tests-factoranalysis.R` algebraic tests for the two-factor analysis model

