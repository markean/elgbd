
<!-- README.md is generated from README.Rmd. Please edit that file -->

# elgbd - Empirical Likelihood for General Block Designs

<!-- badges: start -->

[![Project Status: Suspended – Initial development has started, but
there has not yet been a stable, usable release; work has been stopped
for the time being but the author(s) intend on resuming
work.](https://www.repostatus.org/badges/latest/suspended.svg)](https://www.repostatus.org/#suspended)
[![R-CMD-check](https://github.com/markean/elgbd/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/markean/elgbd/actions/workflows/check-standard.yaml)
<!-- badges: end -->

## Overview

The R package **elgbd** performs hypothesis testing for general block
designs with empirical likelihood. The core computational routines are
implemented with the
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) C++
library and [RcppEigen](https://cran.r-project.org/package=RcppEigen)
interface, with OpenMP for parallel computation. Details of the testing
procedures are given in [Kim, MacEachern, and Peruggia
(2021)](https://arxiv.org/abs/2112.09206).

## Installation

You can install the latest development version from
[Github](https://github.com/markean/elgbd).

``` r
# install.packages("devtools")
devtools::install_github("markean/elgbd")
```

## Usage

``` r
library(elgbd)
# analysis of variance
el_aov(formula = Sepal.Length ~ Species, data = iris)
#> Call:
#> el_aov(formula = Sepal.Length ~ Species, data = iris)
#> 
#> minimizer:
#> 5.9595 5.9595 5.9595
#> 
#> statistic:
#> 94.9054
```
