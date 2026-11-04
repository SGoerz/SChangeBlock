# SChangeBlock

## Spatial Structural Change Detection by an Analysis of Variability Between Blocks of Observations
<!-- badges: start -->
[![R-CMD-check](https://github.com/SGoerz/SChangeBlock/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SGoerz/SChangeBlock/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/SChangeBlock)](https://CRAN.R-project.org/package=SChangeBlock)
<!-- badges: end -->

## Installation

The official version of SChangeBlock can be installed using the R Console:

``` r
install.packages("SChangeBlock")
```

You can install the development version of glmSTARMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SGoerz/SChangeBlock")
```

This `R` package proves methods to detect structural changes in time series or random fields (spatial data). Focus is on the detection of abrupt changes or trends in independent data, but the package also provides a function to de-correlate data with dependence. The functions are based on the test suggested in Schmidt (2024) ([doi:10.3150/23-BEJ1686](https://doi.org/10.3150%2F23-BEJ1686)) and the work in Görz and Fried (2026) ([doi:10.48550/arXiv.2512.11599](https://doi.org/10.48550%2FarXiv.2512.11599)).
