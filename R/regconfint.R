#' @title Estimates bootstrap confidence interval for regularized regression
#' @description \code{regconfint} function estimates confidence intervals for regularized regression coefficients through bootstrapping. The function supports various regression families. For elastic-net, the function supports automatic alpha selection when multiple alpha values are provided, and it also allows for parallel computation to optimize performance.
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @param dataset A data frame containing the data used for the estimation of confidence interval. It should include both the explanatory and target variables specified in expv and tarv.
#' @param expv A character vector specifying the names of the explanatory variables (predictors) in the dataset.
#' @param tarv A character string specifying the name of the target variable (response variable) in the dataset.
#' @param itr An integer specifying the number of bootstrap iterations to perform.
#' @param al A numeric string or vector of alpha values for regularized regression: 0 corresponds to ridge regression; 1 corresponds to lasso regression; Values between 0 and 1 corresponds to elastic-net regression.
#' @param sed A numeric value to set the random seed for reproducibility.
#' @param fam A character string specifying the regression family. Supported families include "gausian", "binomial", "poisson", "multinomial", "cox", and "mgaussian." For more information, see R documentation of glmnet.
#' @param lin A character strig specifying the link function for the specified family.
#' @param typ A character string specifying the type of confidence interval to compute. Default is "all". Supported types are those provided by the boot.ci function, such as "norm", "basic", "perc", and "bca". For more information, see the R documentation of boot.
#' @param paral A character string ("Y" or "N") indicating whether to use parallel computation when multiple alpha values are provided for elastic-net. Default is "N".
#' @return return the bootstrap confidence interval from glmnet object
#' @export
#' @examples
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
<<<<<<< HEAD
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
#' # Percentile bootstrap confidence intervals for ridge
#' for (i in 1:num_var) { bt <- boot.ci(ridge_RR,type="perc",index=i)
#' print(c("CI of variable_",i))
#' print(c(bt$t0,bt$percent[4:5])) }
#'
#' #
=======
#' # Ridge regression
#' ridge_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=0, sed=1234, fam="gaussian", lin="identity")
#' print(ridge_RD)
#' # Install boot package if necessary
#' bt_1 <- boot.ci(ridge_RD,type="perc",index=1)
#' print(c("CI of variable_",1))
#' print(c(bt_1$t0,bt_1$percent[4:5]))
#'
#' # Lasso regression
#' # lasso_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=1, sed=1234, fam="gaussian", lin="identity")
#' # Elastic-net
#' # El_RD <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=seq(0.01,0.99,0.01), sed=1234, fam="gaussian",lin="identity", paral="N")
#' # Percentile bootstrap confidence intervals for ridge
#'
#' #Poisson with log-link
#' # Ridge regression
#' ridge_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=0, sed=1234, fam="poisson", lin="log")
#' print(ridge_RR)
#' # Install boot package if necessary
#' bt_2 <- boot.ci(ridge_RR,type="perc",index=1)
#' print(c("CI of variable_",1))
#' print(c(bt_2$t0,bt_2$percent[4:5]))
#'
#' # Lasso regression
#' # lasso_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=1, sed=1234, fam="poisson", lin="log")
#' # Elastic-net
#' # El_RR <- regconfint(dataset, expv=c("x1", "x2", "x3", "x4", "x5", "x6", "x7"), tarv="y", itr=100, al=seq(0.01,0.99,0.01), sed=1234, fam="poisson",lin="log", paral="N")
#' # Percentile bootstrap confidence intervals for ridge
>>>>>>> 74e252fe13eb00abf1517062a5862e7948ca12d2
#'
#' # If paral = "Y", parallel computing will be implemented based on the number of cores available in your computer

regconfint <- function(dataset,expv,tarv,itr,al,sed,fam,lin,typ="all",paral="N")
{if (length(al) == 1){
  get_RR <- function(d,i){
    x1 <- dataset[i,] %>%
      dplyr::select(all_of(expv))
    x <- as.matrix(x1)
    y1 <- dataset[i,] %>%
      dplyr::select(all_of(tarv))
    y <- unlist(y1)

    family_obj <- switch(fam,
                         "gaussian" = gaussian(link = lin),
                         "binomial" = binomial(link = lin),
                         "poisson" = poisson(link = lin),
                         "multinomial" = multinomial(link = lin),
                         "cox" = cox(link = lin),
                         "mgaussian" = mgaussian(link = lin),
                         stop("Invalid fam argument")
    )

    cvfit_1 <- glmnet::cv.glmnet(x, y, family = family_obj,alpha = al)
    s <- cvfit_1$lambda.min
    fit <- glmnet::glmnet(x, y, family = family_obj,alpha = al, lambda = s)
    r <- exp(fit$beta)
    return(r[1:length(expv)])
  }
  set.seed(sed)
  boot_out <- boot::boot(dataset, statistic = get_RR, R = itr)
  # result <- boot::boot.ci(boot_out,type=typ)
  # result
  boot_out
} else if(length(al) > 1 & paral=="Y")
  {get_alpha <- function(){
  mse.df <- NULL
  num_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(num_cores - 2)  # 使用するコア数に合わせて変更
  doParallel::registerDoParallel(cl)
  foreach::foreach(i = 1:length(alpha), .combine = rbind, .packages = 'glmnet',.export = c("fam", "lin","x","y","al") ) %dopar% {
    family_obj <- switch(fam,
                         "gaussian" = gaussian(link = lin),
                         "binomial" = binomial(link = lin),
                         "poisson" = poisson(link = lin),
                         "multinomial" = multinomial(link = lin),
                         "cox" = cox(link = lin),
                         "mgaussian" = mgaussian(link = lin),
                         stop("Invalid fam argument")
    )
    m <- glmnet::cv.glmnet(x,y,family = family_obj,alpha = al[i])
    data.frame(alpha = al[i],mse = min(m$cvm))
  }
  parallel::stopCluster(cl)
  best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
  return(best.alpha)
}
get_RR <- function(d,i){
  x1 <- dataset[i,] %>%
    dplyr::select(all_of(expv))
  x <- as.matrix(x1)
  y1 <- dataset[i,] %>%
    dplyr::select(all_of(tarv))
  y <- unlist(y1)

  family_obj <- switch(fam,
                       "gaussian" = gaussian(link = lin),
                       "binomial" = binomial(link = lin),
                       "poisson" = poisson(link = lin),
                       "multinomial" = multinomial(link = lin),
                       "cox" = cox(link = lin),
                       "mgaussian" = mgaussian(link = lin),
                       stop("Invalid fam argument")
  )

  best_alpha <- get_alpha()
  cvfit <- glmnet::cv.glmnet(x, y, family = family_obj,alpha = best_alpha)
  s <- cvfit$lambda.min
  fit <- glmnet::glmnet(x, y, family = family_obj,alpha = best_alpha, lambda = s)
  r <- exp(fit$beta)
  return(r[1:length(expv)])
}
set.seed(sed)
boot_out <- boot::boot(dataset, statistic = get_RR, R = itr)
# result <- boot::boot.ci(boot_out,type=typ)
# result
boot_out
} else {
  get_alpha <- function(){
    mse.df <- NULL
    for (i in 1:length(al)) {
      family_obj <- switch(fam,
                           "gaussian" = gaussian(link = lin),
                           "binomial" = binomial(link = lin),
                           "poisson" = poisson(link = lin),
                           "multinomial" = multinomial(link = lin),
                           "cox" = cox(link = lin),
                           "mgaussian" = mgaussian(link = lin),
                           stop("Invalid fam argument")
      )
      m <- glmnet::cv.glmnet(x,y,family = family_obj,alpha = al[i])
      mse.df <- rbind(mse.df, data.frame(alpha = al[i],
                                         mse = min(m$cvm)))
    }
    best.alpha <- mse.df$alpha[mse.df$mse == min(mse.df$mse)]
    return(best.alpha)
  }
  get_RR <- function(d,i){
    x1 <- dataset[i,] %>%
      dplyr::select(all_of(expv))
    x <- as.matrix(x1)
    y1 <- dataset[i,] %>%
      dplyr::select(all_of(tarv))
    y <- unlist(y1)

    family_obj <- switch(fam,
                         "gaussian" = gaussian(link = lin),
                         "binomial" = binomial(link = lin),
                         "poisson" = poisson(link = lin),
                         "multinomial" = multinomial(link = lin),
                         "cox" = cox(link = lin),
                         "mgaussian" = mgaussian(link = lin),
                         stop("Invalid fam argument")
    )

    best_alpha <- get_alpha()
    cvfit <- glmnet::cv.glmnet(x, y, family = family_obj,alpha = best_alpha)
    s <- cvfit$lambda.min
    fit <- glmnet::glmnet(x, y, family = family_obj,alpha = best_alpha, lambda = s)
    r <- exp(fit$beta)
    return(r[1:length(expv)])
  }
  set.seed(sed)
  boot_out <- boot::boot(dataset, statistic = get_RR, R = itr)
  # result <- boot::boot.ci(boot_out,type=typ)
  # result
  boot_out
}
}

