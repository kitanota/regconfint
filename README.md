
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

#Data set creation
 n <- 500
 x1 <- rnorm(n, mean = 50, sd = 10)
 x2 <- rnorm(n, mean = 100, sd = 15)
 x3 <- rnorm(n, mean = 30, sd = 5)
 x4 <- rnorm(n, mean = 70, sd = 20)
 x5 <- rbinom(n, 1, 0.5)
 x6 <- rbinom(n, 1, 0.3)
 x7 <- rbinom(n, 1, 0.7)
 x <- data.frame(x1, x2, x3, x4, x5, x6, x7)
 num_var <- ncol(x)

 #Binary outcome
 y_1 <- 3* x1 - 2 * x2 + x3 + 0.5 * x4 + x5 * 3 + x6 * (-4) + x7 * 6 + rnorm(n, mean = 0, sd = 20)
 prob <- 1 / (1 + exp(-y_1))
 prob[prob>1] <- 1
 prob[prob<0] <- 0
 y <- rbinom(n, 1, prob)
 dataset <- data.frame(x,y)

 #gaussian with identity-link
 # Ridge regression
 ridge_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=0, sed=1234, fam="gaussian", lin="identity")
 # Lasso regression
 lasso_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=1, sed=1234, fam="gaussian", lin="identity")
 # Elastic-net
 El_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=seq(0.01,0.99,0.01), sed=1234, fam="gaussian",lin="identity", paral="N")
 # Percentile bootstrap confidence intervals for ridge
 for (i in 1:num_var) {
 bt <- boot.ci(ridge_RD,type="perc",index=i)
 print(c(Pena_name, "CI of variable_",i))
 print(c(bt$t0,bt$percent[4:5]))
 }

 #Poisson with log-link
 # Ridge regression
 ridge_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=0, sed=1234, fam="poisson", lin="log")
 # Lasso regression
 lasso_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=1, sed=1234, fam="poisson", lin="log")
 # Elastic-net
 El_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=200, al=seq(0.01,0.99,0.01), sed=1234, fam="poisson",lin="log", paral="N")
 # Percentile bootstrap confidence intervals for ridge
 for (i in 1:num_var) {
 bt <- boot.ci(ridge_RR,type="perc",index=i)
 print(c("CI of variable_",i))
 print(c(bt$t0,bt$percent[4:5]))
 }
# If paral = "Y", parallel computing will be implemented based on the number of cores available in your computer

``` r
library(regconfint)
## basic example code
```

