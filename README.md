
# regconfint

<!-- badges: start -->
<!-- badges: end -->

The goal of regconfint is to calculate confidence intervals for regularized regression coefficients through bootstrapping. The function supports various regression families. For elastic-net, the function supports automatic alpha selection when multiple alpha values are provided, and it also allows for parallel computation to optimize performance.

## Installation

You can install the development version of regconfint from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kitanota/regconfint")
```

## Example

#Ridge regression
#regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=0, sed=1, fam=c("gaussian",
#"binomial", "poisson", "multinomial", "cox", "mgaussian")poisson", lin="log")

#Lasso regression
#regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=1, sed=1, fam="poisson", lin="log")

#Elastic-net with fixed alpha
#regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=0.5, sed=1, fam="poisson", lin="log")

#Elastic-net with sequence of alpha
#regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=seq(0.01,0.99,0.01), sed=1, fam="poisson", lin="log", paral="N")

``` r
library(regconfint)
## basic example code
```

