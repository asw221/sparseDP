
# An Introduction to the 'sparseDP' package

### Installation

To install the development version of the `sparseDP` package from GitHub,
we recommend running the following in R (version 3.0 or higher),

```R
if (!suppressWarnings (library(devtools, logical = TRUE)))
  install.packages ("devtools")
library (devtools)
install_github("asw221/sparseDP", dependencies = TRUE)
```


If the package does not install out of the box (especially on OSX
operating systems), your compiler may not support OpenMP. Try cloning
the repository and editing the `Makevars` file in `src` as instructed
in the comments in that file.

### Modeling framework
`sparseDP` is intended to provide a solution the homogeneity pursuit
problem from a Bayesian regression perspective, by placing a
horseshoe-centered Dirichlet process prior on the regression
coefficients. The basic idea is to sample the coefficients from an
approximately sparse discrete mixture, promoting coefficient sparse
clustering. For a general regression problem, the prior hierarchy is
laid out below in a stick breaking process representation.


<p align="center"><img src="https://latex.codecogs.com/svg.latex?\large&space;\beta_j&space;|&space;\omega,&space;\theta^*&space;\sim&space;\sum_{h&space;=&space;1}^\infty&space;\omega_h&space;\delta_{\theta^*_h}(\beta_j)" title="\large \beta_j | \omega, \theta^* \sim \sum_{h = 1}^\infty \omega_h \delta_{\theta^*_h}(\beta_j)" /></p>

<p align="center"><img
src="https://latex.codecogs.com/svg.latex?\large&space;\theta^*_h&space;|&space;\tau^{2},&space;\lambda_h^{2}&space;\sim&space;\mathcal{N}(0,&space;\tau^{2}&space;\lambda_h^{2})"
title="\large \theta^*_h | \tau^{2}, \lambda_h^{2} \sim \mathcal{N}(0,
\tau^{2} \lambda_h^{2})" /></p>

<p align="center"><img src="https://latex.codecogs.com/svg.latex?\large&space;\omega&space;\sim&space;\text{stick}(\alpha),&space;\quad&space;\lambda_h&space;\sim&space;C^&plus;(0,&space;1),&space;\quad&space;\tau&space;\sim&space;\pi(\tau)" title="\large \omega \sim \text{stick}(\alpha), \quad \lambda_h \sim C^+(0, 1), \quad \tau \sim \pi(\tau)" /></p>


Where &delta;<sub>z</sub>(&middot;) is the Dirac measure centered on
z, and C<sup>+</sup>(&mu;, &gamma;) denotes the half-Cauchy density
with location &mu; and scale &gamma;.



