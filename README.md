## TestGLTM
R package for testing the goodness-of-fit of Gaussian latent tree models. 
- More information via ?TestGLTM.
- Helper functions are written in `Rcpp`.
- Algebraic tests are implemented in the file `R/tests.R` and the likelihood-ratio test is implemented in `R/LR_EM.R`.
Install via `install.packages(TestGLTM_1.0.tar.gz, repos=NULL, type="source")`.

## simulations
R-files to compute empirical test sizes and empirical powers via Monte Carlo simulations. 
- Applies the tests in the package TestGLTM.
- Simulations are parallelized to run on a cluster using `foreach`.

## results
The results of the simulations. For each experiment there is a rds-file containing the data of the experiment and a pdf-file containing the corresponding plot. The directory results/final_plots contains the plots displayed in the thesis.

## drawings
Drawings of trees generated via draw.io. The drawings are displayed in the thesis.
