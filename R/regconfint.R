#' @title Estimate confidence interval for regularized regression
#' @description \code{regconfint} regconfint can estimate confidence interval for regularized regression
#' @importFrom dplyr select
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @param dataset data.frame for estimation of confidence intervals
#' @param expv vector of explanatory variables
#' @param tarv target variable
#' @param itr number of iteration
#' @param al alpha for regularized regression:0(Ridge),1(Lasso) or sequence of alpha for elastic-net
#' @param sed seed for bootstrap
#' @param fam Either a character string representing one of the built-in families, or else a glm() family object. For more information, see R documentation of glmnet.
#' @param lin link function of the model
#' @param typ type of intervals required. Default is "all", For more information, see the R documentation of boot.
#' @return return the bootstrap confidence interval from glmnet object
#' @export
#' @examples
#' #Ridge regression
#' #regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=0, seed=1, fam=c("gaussian",
#' #"binomial", "poisson", "multinomial", "cox", "mgaussian")poisson", lin="log")
#'
#' #Lasso regression
#' #regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=1, seed=1, fam="poisson", lin="log")
#'
#' #Elastic-net
#' #regconfint(dataset, expv=c("x1","x2"), tarv="y", itr=50, al=seq(0.01,0.99,0.01), seed=1, fam="poisson",
#' # lin="log")



regconfint <- function(dataset,expv,tarv,itr,al,sed,fam,lin,typ="all")
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
  result <- boot::boot.ci(boot_out,type=typ)
  result
}else{
  get_alpha <- function(){
    mse.df <- NULL
    for (i in 1:length(al)) {
      m <- glmnet::cv.glmnet(x,y,family = fam(link=lin),alpha = al[i])
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
  result <- boot::boot.ci(boot_out,type=typ)
  result
}
}
