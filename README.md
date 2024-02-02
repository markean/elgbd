
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
