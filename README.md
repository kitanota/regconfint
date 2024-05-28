
# regconfint

<!-- badges: start -->
<!-- badges: end -->

The goal of regconfint is to ...

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
#regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=seq(0.01,0.99,0.01), sed=1, fam="poisson",
# lin="log", paral="N")

# If paral = "Y", parallel computing will be implemented based on the number of cores available
# in your computer

``` r
library(regconfint)
## basic example code
```

