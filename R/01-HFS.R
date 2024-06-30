#' hybridFS: Hybrid Feature Screening for high-dimensional data.
#'
#' @param X An n by p numerical matrix. Each row corresponds to an observation.
#' @param y The response numerical vector.
#' @param tau A numerical cutoff between 0 and 1 for HFS. Default 0.5.
#' @param prop The maximal proportion of kept features. Default 0.9.
#' @param isAll If set to \code{TRUE}, set to a stricter standard in screening. Default FALSE.
#' @param utilities Default \code{list()}. If not specified, HFS will compute KPC
#'  and R-squared from polynomial regression instead.
#' @param M The order of the polynomial regression if \code{utilities} not specified. Default 3.
#'
#'
#' @return Returns the index of the selected variables and their variable scores.
#'
#' @importFrom KPC KPCRKHS
#' @importFrom isotree isolation.forest
#' @importFrom dHISC dhsic
#' @importFrom parallel makeCluster
#'
#' @export
hybrid.corr <- function(X, y, tau=0.5, prop=0.9, isAll=FALSE, utilities=list(), M=3){
  n <- dim(X)[1]
  if(n != length(y)){
    stop("The length of y does not match the row number of X")
  }
  p <- dim(X)[2]
  v <- length(utilities)
  type <- ifelse(length(unique(y))==2, "binary-classification", "regression")
  
  # use parallel computing to save time
  print(parallel::detectCores())
  if(parallel::detectCores() > 1){
    cl <- parallel::makeCluster(min(parallel::detectCores()-1, 8))
  }else{
    cl <- parallel::makeCluster(parallel::detectCores())
  }
  doParallel::registerDoParallel(cl)
  # compute the default utilities if not specified
  if(v == 0){
    v = 2
    corrs <- matrix(NA, p, 2)
    corrs[, 1] <- foreach::foreach(i=1:p, .combine="c") %dopar% {
      x <- X[, i]
      x_powers <- matrix(nrow=length(y), ncol=M)
      for(i in seq(M)){
        x_powers[, i] <- x ^ i
      }
      if (type == "binary-classification"){
        model <- glm(y ~ x_powers, family = "binomial")
        null_model <- glm(y ~ 1, family = "binomial")
        PseudoR2 <- 1 - (logLik(model) / logLik(null_model))
        return(unname(PseudoR2))
      }else {
        return(summary(lm(y~x_powers))$adj.r.squared)
      }
    }
    corrs[, 2] <- foreach::foreach(i=1:p, .combine="c") %dopar% {
      rho <- NaN
      try(rho <- KPC::KPCRKHS(Y=y, Z=X[, i]))
    }
    if(sum(is.na(corrs[, 2])) >= floor(p*0.1)){
      corrs[, 2] <- foreach::foreach(i=1:p, .combine="c") %dopar% dHSIC::dhsic(X=X[, i], Y=y)$dHSIC
    }
  }else{
    # otherwise compute the user specified utility functions
    corrs <- matrix(NA, p, v)
    for(i in seq(v)){
      u_fun <- utilities[[i]]
      corrs[, i] <- sapply(1:p, function(i) u_fun(X[, i], y))
    }
  }
  suppressWarnings(parallel::stopCluster(cl))
  
  # compute the feature-wise isolation score
  fit.isotree <- isotree::isolation.forest(corrs, ndim=v, ntrees=100, nthreads=1)
  score <- predict(fit.isotree, corrs)
  cutoff <- apply(corrs, 2, function(x) quantile(x, probs=prop, na.rm=TRUE))
  if (isAll){
    cond <- apply(corrs, 1, function(x) all(x >= cutoff))
  }else{
    cond <- apply(corrs, 1, function(x) TRUE %in% c(x>=cutoff))
  }
  idx <- which(cond == TRUE & score >= tau)
  return(list(idx=idx, corrs=corrs, score=score, type=type))
}

