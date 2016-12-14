## Wedge probability computation

This `R` package is devoted to the computation of wedge probabilities. It includes two efficient functions using pure `R` vectorization and `Rcpp` interfacing `C++` code directly in the `R` system.

### Installation

Currently, the installation of the package is made using the `devtools` package:

```{bash}
devtools::install_github("rcqls/wedge")
```

### Interested in faster computation of wedge probabilities?

A companion `R` package [`wedgeParallel`](http://github.com/rcqls/wedgeParallel) is also available. It takes advantage of the multicore inside our modern computers. To install it, simply execute the following line in a terminal:

```{bash}
devtools::install_github("rcqls/wedgeParallel")
```
