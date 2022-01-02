#!/bin/bash

# Run file via ". dev-Rsetup.sh". Then the current shell context will be used.

# Set paths to R libraries
R_LIBS_USER="/dss/dsshome1/lxc0D/ge73wex3/R/debian10/3.6/"
export R_LIBS_USER
R_LIBS_SITE="/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/haswell/r/3.6.3-gcc-pv52rxc/rlib/R/library"
export R_LIBS_SITE

# Load R and compatible gcc compiler and gsl
module load r/3.6.3-gcc8-mkl
module unload intel-mpi
module unload intel
module load gcc/8.4.0
module load gsl/2.5-gcc8
