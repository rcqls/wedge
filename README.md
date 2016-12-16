# Wedge probability computation

This `R` package `wedge` is devoted to the computation of wedge probabilities. It includes two functions using pure `R` vectorization and `Rcpp` (interfacing `C++` code directly in the `R` system).

## Installation

Installation `wedge` using the `devtools` package:

```{bash}
devtools::install_github("rcqls/wedge")
```

Binaries for Windows and MacOS are available upon request

## Interested in parallel computation?

A companion `R` package [`wedgeParallel`](http://github.com/rcqls/wedgeParallel) is also available. It takes full advantage of the multicore architecture of your computer.
To install it, run the following line in a terminal:

```{bash}
devtools::install_github("rcqls/wedgeParallel")
```
