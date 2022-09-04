## TestGGM
R package for testing the goodness-of-fit of Gaussian graphical models. Supports two types of models:
- Gaussian latent tree models
- Two-factor analysis model

Install via `devtools::install_github("NilsSturma/TestGGM")`. The package implements the test using incomplete U-statistics described in 

Sturma, N., Drton, M., and Leung, D. (2022).
Testing Many and Possibly Singular Polynomial Constraints.
[arXiv](https://arxiv.org/abs/2208.11756).

Moreover, it supports similar testing strategies by applying different schemes of forming averages over the kernel estimators.

Most important files:
- `R/tests.R`: algebraic tests for (arbitrary) Gaussian latent tree models 
- `R/LR_EM.R`: likelihood-ratio test for (arbitrary) Gaussian latent tree models using the EM algorithm
- `R/tests-factoranalysis.R` algebraic tests for the two-factor analysis model

