# MESuSiE

## Multiple Ancestry Sum of Single Effect Model (MESuSiE)

![](MESuSiE_Overview.png)

MESuSiE relies on GWAS summary statistics from multiple ancestries, properly accounts for the LD structure of the local genomic region in multiple ancestries, and explicitly models both shared and ancestry-specific causal signals to accommodate causal effect size similarity as well as heterogeneity across ancestries. MESuSiE outputs posterior inclusion probability of variant being shared or ancestry-specific causal variants. 
MESuSiE is implemented as an open-source R package, freely available at [Zhou lab](https://www.xzlab.org/software.html).

## Installation

To install the latest version of the MESuSiE package from GitHub, run
the following code in R:

```R
library(devtools)
install_github("borangao/MESuSiE")
```

This command should automatically install all required packages if
they are not installed already.

## Quick Start

See [Tutorial](https://borangao.github.io/meSuSie_Analysis/) for detailed documentation and examples.

## Reproduce

### Optional covariance EM

The author covariance optimizer remains the default (`optim_method="optim"`).
An optional `optim_method="em"` uses permutation-equivariant covariance EM
updates and stable Gaussian likelihood/posterior arithmetic:

```r
fit <- meSuSie_core(LD_list, summ_stat_list, L=3, optim_method="em",
                   em_max_iter=100, em_tol=1e-9)
fit$converged       # Outer stopping only
fit$em_maxiter      # Inner calls that exhausted the iteration budget
fit$em_status      # Last inner status for each effect
```

Both optimizers now require a finite `0 <= ELBO change < tol` for outer
convergence, following [updated SuSiE](https://github.com/stephenslab/susieR/commit/b5c71c1e7c11b33f82ec230fbaa0c386fd90c55b).
A decrease is not convergence. `converged=FALSE` distinguishes a returned fit
that reaches `max_iter`. The default `tol` is unchanged at .001.
Finite-budget EM is not identical to the original optimizer; outer convergence
does not imply converged inner maximization, a global optimum, or corrected LD.
The EM update floors covariance eigenvalues at 1e-12 and retains the previous
covariance if the proposed update lowers the evaluated likelihood.

Base-R numerical and population-permutation regressions can be run after
installation with `Rscript tests/numerical-regressions.R`. The optional kernels
reuse MIT-licensed credtools code (notice in `inst/CREDTOOLS_LICENSE`).

See [Repository](https://doi.org/10.5281/zenodo.8411004) for reproducing the simulation and real data analysis in the manuscript. 


## Issues
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R! 

How to cite `MESuSiE`
-------------------
Boran Gao, Xiang Zhou#. MESuSiE enables scalable and powerful multi-ancestry fine-mapping of causal variants in genome-wide association studies.
