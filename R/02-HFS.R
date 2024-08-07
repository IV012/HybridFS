#' hybridFS: Hybrid Feature Screening for high-dimensional data.
#'
#' @param X An n by p numerical matrix. Each row corresponds to an observation.
#' @param y The numerical response vector.
#' @param tau A numerical cutoff between 0 and 1 for HFS. Default 0.5.
#' @param prop The maximal proportion of kept features. Default 0.9.
#' @param isAll If set to \code{TRUE}, set to a stricter standard in screening. Default FALSE.
#' @param utilities Default \code{list()}. If not specified, HFS will compute HSIC
#'  and R-squared from polynomial regression instead.
#' @param M The order of the polynomial regression if \code{utilities} not specified. Default 3.
#'
#'
#' @return Returns list of variable indexes \code{idx} and HFS scores \code{score}.
#'
#' @importFrom isotree isolation.forest
#' @importFrom dHSIC dhsic
#' @importFrom parallel makeCluster
#' @importFrom foreach foreach
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats logLik
#'
#' @export
hybrid.corr <- function(X, y, tau=0.5, prop=0.9, isAll=FALSE, utilities=list(), M=3){
  n <- dim(X)[1]
  if(n != length(y)){
    stop("The length of y does not match the row number of X")
  }
  p <- dim(X)[2]
  if (p < 10){
    stop("Too few features as input. HFS accepts at least 10 features.")
  }
  v <- length(utilities)
  type <- ifelse(length(unique(y))==2, "binary-classification", "regression")

  # use parallel computing to save time
  # print(parallel::detectCores())
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
        model <- stats::glm(y ~ x_powers, family = "binomial")
        null_model <- stats::glm(y ~ 1, family = "binomial")
        PseudoR2 <- 1 - (stats::logLik(model) / stats::logLik(null_model))
        return(unname(PseudoR2))
      }else {
        return(summary(stats::lm(y~x_powers))$adj.r.squared)
      }
    }
    corrs[, 2] <- foreach::foreach(i=1:p, .combine="c") %dopar% dHSIC::dhsic(X=X[, i], Y=y)$dHSIC
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
  if (length(idx) <= 10){
    idx <- which(rank(-score) <= 10)
  }
  return(list(idx=idx, corrs=corrs, score=score, type=type))
}

#' hybridFS: Hybrid Feature Screening for high-dimensional data.
#'
#' @param X An n by p numerical matrix. Each row corresponds to an observation.
#' @param y The numerical response vector.
#' @param object An object from \code{hybrid.corr}.
#' @param val_idx The index of validation data. Default \code{NULL}.
#' @param tau_set Candidate taus in an increasing order.
#'  Default (0.5, 0.55, 0.6, 0.65, 0.7).
#' @param val_ratio If \code{val_idx} is not specified, randomly split the ratio
#'  of data to a validation set. Default 0.2.
#' @param mod_fun A training function passed to \code{permfit}.
#'  Default \code{mod_permfit}.
#' @param predict_fun A prediction function passed to \code{permfit}.
#'  Defau \code{predict_mod_permfit}.
#' @param method Models, including \code{random_forest} for random forests,
#'  \code{lasso} for linear/logistic lasso, \code{linear} for
#'  linear/logistic regression,
#'   \code{svm} for svms with Gaussian kernels.
#' @param alpha The cutoff for p-values. Default 0.1.
#' @param prop The maximal proportion of kept features. Default 0.9.
#' @param isAll If set to \code{TRUE}, set to a stricter standard in screening.
#'  Default FALSE.
#' @param ... Other parameters passed to \code{permfit}.
#'
#' @return The list of variable indexes \code{idx} and best tau \code{tau.best}.
#'
#' @export
tune.tau <- function(X, y, object, val_idx = NULL,
                     tau_set = c(0.5, 0.55, 0.6, 0.65, 0.7), val_ratio = 0.2,
                     mod_fun = mod_permfit, predict_fun = predict_mod_permfit,
                     method = c("random_forest", "lasso",
                                "linear", "svm", "xgboost")[1],
                     alpha = 0.1, prop = 0.95, isAll = FALSE, ...){
  if (length(tau_set) <= 1) {
    stop("The function 'tune.tau' only accept multiple taus. Please use 'hybridFS' instead.")
  }
  if (dim(X)[2] < 10){
    stop("Too few features as input. HFS accepts at least 10 features.")
  }
  tau_set <- unique(tau_set[order(tau_set)])
  # split the data into training and testing sets
  n <- length(y)
  if (is.null(val_idx)){
    val_idx <- sample(n, n*min(val_ratio, 0.5))
  }
  if (length(val_idx) <= 5) {
    warning("Too Few validation Samples to Tune Tau.")
    return(list(idx=object$idx, tau.best=tau_set[1], hfs.object=object, method=method))
  }
  X_t <- X[val_idx, ]
  y_t <- y[val_idx]
  X <- X[-val_idx, ]
  y <- y[-val_idx]
  fit.type <- object$type
  if (fit.type == "binary-classification"){
    y <- as.factor(y)
    y_t <- as.factor(y_t)
  }
  # condition on the quantile
  score <- object$score
  cutoff <- apply(object$corrs, 2, function(x) quantile(x, probs=prop, na.rm=TRUE))
  if (isAll){
    cond <- apply(object$corrs, 1, function(x) all(x >= cutoff))
  }else{
    cond <- apply(object$corrs, 1, function(x) TRUE %in% c(x >= cutoff))
  }
  # fast group permfit to roughly determine the cutoff
  min.tau <- max(tau_set)
  for (tau in tau_set){
    fea_idx <- which(score >= tau & cond == TRUE)
    if (length(fea_idx) <= n / 5) {
      min.tau <- tau
      break
    }
  }
  tau_set <- tau_set[tau_set >= min.tau]
  score <- score[fea_idx]
  train <- importDnnet(x = X[, fea_idx], y = y)
  validate <- importDnnet(x = X_t[, fea_idx], y = y_t)
  # create blocks of indexes within the intervals of tau.set
  pathway <- list()
  for (i in seq(length(tau_set) - 1)){
    block.name <- paste("Block", tau_set[i], sep="")
    block.idx <- which(score >= tau_set[i] & score < tau_set[i+1])
    if (length(block.idx) > 1){
      pathway[[block.name]] <- block.idx
    }
  }
  block.imps <- permfit(
    train = train, validate = validate, k_fold = 0,
    active_var = c(), mod_fun = mod_fun,
    predict_fun = predict_fun, pathway_list = pathway, method = method,
    ...)@block_importance$importance_pval_x
  block.keys <- which(block.imps <= alpha)
  if(length(block.keys) > 0){
    pathway <- pathway[block.keys]
    tau.best <- as.numeric(gsub("Block", "", names(pathway)[1]))
  }else{
    tau.best <- max(tau_set)
  }
  idx <- fea_idx[which(score >= tau.best)]
  if (length(idx) <= 10){
    idx <- which(rank(-score) <= 10)
  }
  return(list(idx=idx, tau.best=tau.best, hfs.object=object, method=method))
}
