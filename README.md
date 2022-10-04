
<!-- README.md is generated from README.Rmd. Please edit that file -->

# elgbd - Empirical Likelihood for General Block Designs

<!-- badges: start -->

[![R-CMD-check](https://github.com/markean/elgbd/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/markean/elgbd/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The R package **elgbd** provides a unified framework for data analysis
with empirical likelihood methods. A collection of functions are
available for basic regression analysis and hypothesis testing. Much of
its functionality and syntax are designed to mimic the corresponding
base R functions. The core routines are written in C++ and utilize
OpenMP for parallelization.

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
