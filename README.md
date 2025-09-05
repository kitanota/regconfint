
<!-- README.md is generated from README.Rmd. Please edit that file -->

# regconfint

<!-- badges: start -->
<!-- badges: end -->

The goal of regconfint is to …

## Installation

You can install the development version of regconfint from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("kitanota/regconfint")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(regconfint)
#' #Data set creation
#' n <- 500
#' x1 <- rnorm(n, mean = 50, sd = 10)
#' x2 <- rnorm(n, mean = 100, sd = 15)
#' x3 <- rnorm(n, mean = 30, sd = 5)
#' x4 <- rnorm(n, mean = 70, sd = 20)
#' x5 <- rbinom(n, 1, 0.5)
#' x6 <- rbinom(n, 1, 0.3)
#' x7 <- rbinom(n, 1, 0.7)
#' x <- data.frame(x1, x2, x3, x4, x5, x6, x7)
#' num_var <- ncol(x)
#'
#' #Binary outcome
#' y_1 <- 3* x1 - 2 * x2 + x3 + 0.5 * x4 + x5 * 3 + x6 * (-4) + x7 * 6 + rnorm(n, mean = 0, sd = 20)
#' prob <- 1 / (1 + exp(-y_1))
#' prob[prob>1] <- 1
#' prob[prob<0] <- 0
#' y <- rbinom(n, 1, prob)
#' dataset <- data.frame(x,y)
#' print(dataset)
#'
#' #gaussian with identity-link
#'
#' # Ridge regression
#' ridge_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=0, sed=1234, fam="gaussian", lin="identity")
#'
#' # Lasso regression
#' lasso_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=1, sed=1234, fam="gaussian", lin="identity")
#'
#' # Elastic-net
#' El_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=seq(0.01,0.99,0.01), sed=1234, fam="gaussian",lin="identity", paral="N")
#'
#' # Percentile bootstrap confidence intervals for ridge
#' for (i in 1:num_var) { bt <- boot.ci(ridge_RD,type="perc",index=i)
#' print(c(Pena_name, "CI of variable_",i))
#' print(c(bt$t0,bt$percent[4:5])) }
#'#'
#' #Poisson with log-link
#'
#' # Ridge regression
#' ridge_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=0, sed=1234, fam="poisson", lin="log")
#'
#' # Lasso regression
#' lasso_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=1, sed=1234, fam="poisson", lin="log")
#' # Elastic-net
#' El_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=seq(0.01,0.99,0.01), sed=1234, fam="poisson",lin="log", paral="N")
#'
#' # If paral = "Y", parallel computing will be implemented based on the number of cores available in your computer
#'
#' # Bias-corrected bootstrap confidence intervals for ridge
#' BC_CI <- bc_ci(ridge_RR,conf=0.95)
#' print(BC_CI)
#'
#' # Percentile bootstrap confidence intervals for ridge
#' # For other types of confidence intervals, modify the type argument. For more details, please see the instruction of boot pakage
#' for (i in 1:num_var) { bt <- boot.ci(ridge_RR,type="perc",index=i)
#' print(c("CI of variable_",i))
#' print(c(bt$t0,bt$percent[4:5])) }
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
