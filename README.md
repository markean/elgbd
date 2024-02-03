
<!-- README.md is generated from README.Rmd. Please edit that file -->

# elgbd

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/markean/elgbd/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/markean/elgbd/actions/workflows/check-standard.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/elgbd)](https://CRAN.R-project.org/package=elgbd)
<!-- badges: end -->

## Overview

elgbd performs hypothesis testing for general block designs with
empirical likelihood. The core computational routines are implemented
with the ‘Eigen’ ‘C++’ library and ‘RcppEigen’ interface, with ‘OpenMP’
for parallel computation. Details of the testing procedures are given in
[Kim, MacEachern, and Peruggia
(2023)](https://doi.org/10.1080/10485252.2023.2206919). This work was
supported by the U.S. National Science Foundation under Grants
No. [SES-1921523](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1921523)
and
[DMS-2015552](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2015552).

## Installation

You can install the latest stable release of elgbd from CRAN.

``` r
install.packages("elgbd")
```

### Development version

You can install the development version of elgbd from GitHub.

``` r
# install.packages("pak")
pak::pak("markean/elgbd")
```

## Usage

``` r
library(elgbd)
set.seed(58234)
# Analysis of variance
el_aov(formula = Sepal.Length ~ Species, data = iris)
#> Call:
#> el_aov(formula = Sepal.Length ~ Species, data = iris)
#> 
#> Minimizer:
#> 5.9595 5.9595 5.9595
#> 
#> Statistic:
#> 94.9054

# All pairwise comparisons
data("clothianidin")
el_pairwise(clo ~ trt | blk, data = clothianidin, B = 1000)
#> 
#>  Empirical Likelihood Multiple Tests
#> 
#> All pairwise comparisons
#> 
#>                   Estimate  Chisq  Lwr.ci  Upr.ci  p.adj    
#> Naked - Fungicide  -1.0525 14.098 -2.1150 -0.2794  0.002 ** 
#> Naked - Low        -1.6794 15.935 -2.6705 -0.6435  0.001 ***
#> Naked - High       -3.1726 42.803 -4.5731 -1.6940 <0.001 ***
#> Fungicide - Low    -0.6269  1.834 -1.6303  0.6448  0.523    
#> Fungicide - High   -2.1201 14.085 -3.2313 -0.7247  0.002 ** 
#> Low - High         -1.4932 11.715 -2.7228 -0.3954  0.003 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> k: 1, level: 0.05, method: AMC, cutoff: 6.6737

# Comparisons with control
el_pairwise(clo ~ trt | blk,
  data = clothianidin, control = "Naked", method = "NB", B = 1000
)
#> 
#>  Empirical Likelihood Multiple Tests
#> 
#> Comparisons with control
#> 
#>                   Estimate Chisq Lwr.ci Upr.ci p.adj   
#> Fungicide - Naked   1.0525 14.10 0.3506 1.9900 0.005 **
#> Low - Naked         1.6794 15.94 0.7586 2.5646 0.004 **
#> High - Naked        3.1726 42.80 1.8489 4.4282 0.003 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> k: 1, level: 0.05, method: NB, cutoff: 6.3291

# Test for equal means
el_test(clo ~ trt | blk, clothianidin,
  lhs = matrix(c(
    1, -1, 0, 0,
    0, 1, -1, 0,
    0, 0, 1, -1
  ), byrow = TRUE, nrow = 3L)
)
#> 
#>  Empirical Likelihood Test
#> 
#> General block designs
#> 
#> Maximum EL estimates:
#> [1] -4.479086 -3.426618 -2.799700 -1.306530
#> 
#> Statistic: 28.055
```
