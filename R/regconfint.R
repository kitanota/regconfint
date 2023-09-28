#' @title Estimating confidence interval for regularized regression
#' @description \code{regconfint} Estimating confidence interval for regularized regression
#' @importFrom dplyr SELECT
#' @importFrom glmnet CV.GLMNET
#' @importFrom glmnet GLMNET
#' @importFrom boot BOOT
#' @importFrom boot BOOT.CI
#' @param dataset dataset for calculating CI
#' @param expv explanatory variables(character)
#' @param tarv target variable(character)
#' @param itr number of iteration
#' @param al alpha for regularized regression:0(Ridge),1(Lasso) or sequence of alpha for elastic-net
#' @param sed seed for bootstrap
#' @param fam distribution of the model
#' @param lin link function of the model
#' @param typ type of intervals required.
#' @return return the bootstrap confidence interval from glmnet object
#' @export
#' @examples
#' #regconfint(dataset=xxxx, expv=c("a","b","c"), tarv="target", itr=1000, al=1, seed=1234, fam="poisson", lin="log")

regconfint <- function(dataset,expv,tarv,itr,al,sed,fam,lin,typ="all")
{if (al == 0|al == 1){
  get_RR <- function(d,i){
    x1 <- dataset[i,] %>%
      dplyr::select(all_of(expv))
    x <- as.matrix(x1)
    y1 <- dataset[i,] %>%
      dplyr::select(all_of(tarv))
    y <- unlist(y1)
    cvfit_1 <- glmnet::cv.glmnet(x, y, family = fam(link=lin),alpha = al)
    s <- cvfit_1$lambda.min
    fit <- glmnet::glmnet(x, y, family = fam(link=lin),alpha = al, lambda = s)
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
    best_alpha <- get_alpha()
    cvfit <- glmnet::cv.glmnet(x, y, family = fam(link=lin),alpha = best_alpha)
    s <- cvfit$lambda.min
    fit <- glmnet::glmnet(x, y, family = fam(link=lin),alpha = best_alpha, lambda = s)
    r <- exp(fit$beta)
    return(r[1:length(expv)])
  }
  set.seed(sed)
  boot_out <- boot::boot(dataset, statistic = get_RR, R = itr)
  result <- boot::boot.ci(boot_out,type=typ)
  result
}
}
