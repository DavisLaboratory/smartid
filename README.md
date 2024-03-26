
# smartid: Scoring and MARker identification method based on modified Tf-IDf

<!-- badges: start -->
[![R-CMD-check](https://github.com/DavisLaboratory/smartid/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/DavisLaboratory/smartid/actions)
[![Coverage status](https://codecov.io/gh/DavisLaboratory/smartid/branch/devel/graph/badge.svg)](https://codecov.io/github/DavisLaboratory/smartid?branch=devel)
[![BioC status](https://bioconductor.org/shields/years-in-bioc/smartid.svg)](https://bioconductor.org/packages/smartid/)
<!-- badges: end -->

An R package for automatically identify group specific signature and score cells based on given gene sets on scRNA data.

smartid is an R package designed for automated identification of signatures of interest. The package is developed for generating specific signature genes from multiple groups based on modified TF-IDF approach, which is good at finding markers for rare populations and pay more attention to the genes without high expression but with high variance across groups.

This package is particularly useful for the marker identification of novel or rare group populations in various biological and medical applications, including cancer research and developmental biology.

[**Check out the standard demonstration.**](https://davislaboratory.github.io/smartid/articles/smartid_Demo.html)

## Installation

You can install the development version of smartid like so:

``` r
# install from GitHub
devtools::install("DavisLaboratory/smartid")
```

smartid can be installed from Bioconductor directly as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("smartid")
```

