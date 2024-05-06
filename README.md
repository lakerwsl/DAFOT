# DAFOT

<!-- Badges begin -->
[![R-CMD-check](https://github.com/lakerwsl/DAFOT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lakerwsl/DAFOT/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/lakerwsl/DAFOT/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/lakerwsl/DAFOT/actions/workflows/test-coverage.yaml)
[![Published in Biometrika](https://img.shields.io/badge/Published%20in-Biometrika-blue)](https://doi.org/10.1093/biomet/asaa061)
[![Published in Biometrika](https://img.shields.io/badge/Published%20in-Biometrika-blue)](https://doi.org/10.1093/biomet/asad075)
<!-- Badges end -->

This is an R package implementing the DAFOT test proposed in Wang, S., Cai, T.T., Li, H. (2020) and the IndDAFOT and ConIndDAFOT tests proposed in Wang, S., Yuan, B., Cai, T.T., Li, H. (2023). The DAFOT test is designed to test the association between the phylogenetic composition of a microbial community and a continuous outcome. The IndDAFOT and ConIndDAFOT tests are designed to test the association between the phylogenetic composition of a microbial community and a continuous outcome with adjustment for covariates.

### Installation

You can install the released version of DAFOT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DAFOT")
```

You can install the development version of DAFOT from [GitHub](https://github.com/lakerwsl/DAFOT) with:

``` r
# install.packages("devtools")
devtools::install_github("lakerwsl/DAFOT")
```

### Quick Start

You can run the core functions on our test data sets to get a sense for how they work. Starting with `DAFOT`:

``` r
library(DAFOT)

test_DAFOT

dafot <- DAFOT(test_DAFOT$P, test_DAFOT$Q, test_DAFOT$tree)

dafot$Stat
dafot$P
dafot$Thre
dafot$Active
```

You can run test data through `IndDAFOT` and `ConIndDAFOT` very similarly:

``` r
ind_dafot <- IndDAFOT(test_IndDAFOT$P, test_IndDAFOT$Y, test_IndDAFOT$tree)
con_ind_dafot <- ConIndDAFOT(test_ConIndDAFOT$P, test_ConIndDAFOT$Y, test_ConIndDAFOT$X, test_ConIndDAFOT$tree, condgen = test_ConIndDAFOT$condgen)
```

Reference on all functions' inputs and outputs can be found in the package documentation.

### References

For `DAFOT` and `SCalculation` functions, please cite:

Shulei Wang, T Tony Cai, Hongzhe Li, Hypothesis testing for phylogenetic composition: a minimum-cost flow perspective, Biometrika, Volume 108, Issue 1, March 2021, Pages 17–36, https://doi.org/10.1093/biomet/asaa061

For `IndDAFOT` and `ConIndDAFOT` functions, please cite:

Shulei Wang, Bo Yuan, T Tony Cai, Hongzhe Li, Phylogenetic association analysis with conditional rank correlation, Biometrika, 2023;, asad075, https://doi.org/10.1093/biomet/asad075