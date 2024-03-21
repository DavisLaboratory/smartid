
# smartid: Scoring and MARker identification method based on modified Tf-IDf

<!-- badges: start -->
[![R-CMD-check](https://github.com/DavisLaboratory/smartid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DavisLaboratory/smartid/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An R package for automatically identify group specific signature and score cells based on given gene sets on scRNA data.

smartid is an R package designed for automated identification of signatures of interest. The package is developed for generating speicifc signature genes from multiple groups based on modified TF-IDF approach, which is good at finding markers for rare populations and pay more attention to the genes without high expression but with high variance across groups.

This package is particularly useful for the marker identification of novel or rare group populations in various biological and medical applications, including cancer research and developmental biology.

## Installation

You can install the development version of smartid like so:

``` r
# install from GitHub
devtools::install("DavisLaboratory/smartid")
```


