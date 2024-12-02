% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regconfint.R
\name{regconfint}
\alias{regconfint}
\title{Estimates bootstrap confidence interval for regularized regression}
\usage{
regconfint(
  dataset,
  expv,
  tarv,
  itr,
  al,
  sed,
  fam,
  lin,
  typ = "all",
  paral = "N"
)
}
\arguments{
\item{dataset}{A data frame containing the data used for the estimation of confidence interval. It should include both the explanatory and target variables specified in expv and tarv.}

\item{expv}{A character vector specifying the names of the explanatory variables (predictors) in the dataset.}

\item{tarv}{A character string specifying the name of the target variable (response variable) in the dataset.}

\item{itr}{An integer specifying the number of bootstrap iterations to perform.}

\item{al}{A numeric string or vector of alpha values for regularized regression: 0 corresponds to ridge regression; 1 corresponds to lasso regression; Values between 0 and 1 corresponds to elastic-net regression.}

\item{sed}{A numeric value to set the random seed for reproducibility.}

\item{fam}{A character string specifying the regression family. Supported families include "gausian", "binomial", "poisson", "multinomial", "cox", and "mgaussian." For more information, see R documentation of glmnet.}

\item{lin}{A character strig specifying the link function for the specified family.}

\item{typ}{A character string specifying the type of confidence interval to compute. Default is "all". Supported types are those provided by the boot.ci function, such as "norm", "basic", "perc", and "bca". For more information, see the R documentation of boot.}

\item{paral}{A character string ("Y" or "N") indicating whether to use parallel computation when multiple alpha values are provided for elastic-net. Default is "N".}
}
\value{
return the bootstrap confidence interval from glmnet object
}
\description{
\code{regconfint} function estimates confidence intervals for regularized regression coefficients through bootstrapping. The function supports various regression families. For elastic-net, the function supports automatic alpha selection when multiple alpha values are provided, and it also allows for parallel computation to optimize performance.
}
\examples{
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

# If you want to reproduce the authors' study, set seed as "1234" with CMV and NCDS dataset
}
