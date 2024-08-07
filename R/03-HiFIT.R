#' HiFIT: A high-dimensional feature importance test.
#'
#' @param X An n by p numerical matrix. Each row corresponds to an observation.
#' @param y The numerical response vector.
#' @param k_fold K-fold cross-fitting. If not, set k_fold to zero.
#' @param n_perm Number of permutations repeated. Default 100.
#' @param mod_fun Function to train the model. Default \code{mod_permfit}.
#' @param predict_fun Function to predict the response.
#'  Default \code{predict_mod_permfit}.
#' @param pathway_list A list of pathways to be jointly tested.
#' @param active_var Index of active variables.
#' @param method Models, including \code{ensemble_dnnet} for ensemble deep
#'   neural networks, \code{random_forest} for random forests, \code{lasso}
#'   for linear/logistic lasso, \code{linear} for linear/logistic regression,
#'   \code{svm} for svms with Gaussian kernels, and \code{dnnet} with single deep
#'   neural network.
#' @param shuffle If shuffle is null, the data will be shuffled for
#'   cross-fitting; if random shuffle is not desired, please provide
#'   a bector of numbers for cross-fitting indices
#' @param utilities Utility functions passed to \code{hybrid.corr}.
#'  Default \code{list()}.
#' @param M The order of polynomial regression used by HFS. Default 3.
#' @param alpha The cutoff of p-value for feature screening. Default 0.1.
#' @param prop The max proportion of screened features. Default 0.9.
#' @param isAll Adopt stricter screening standard. Default FALSE.
#' @param ... Additional parameters passed to each method.
#'
#' @return Returns a PermFIT object.
#'
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats var
#'
#' @export

hifit <- function(X, y, k_fold = 5,
                    n_perm = 100, mod_fun = mod_permfit,
                    predict_fun = predict_mod_permfit,
                    pathway_list = list(),
                    active_var = NULL,
                    method = c("random_forest", "lasso",
                               "linear", "svm", "xgboost")[1],
                    shuffle = NULL, utilities = list(),
                    M = 3, alpha = 0.1, prop = 0.9,
                    isAll=FALSE, ...) {
  X <- as.matrix(X)
  y <- as.vector(y)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (n != length(y)){
    stop("The number of rows of X should be equal to the length of y.")
  }
  n_pathway <- ifelse(length(pathway_list) > 0, length(pathway_list), 1)
  p_score <- array(NA, dim = c(n_perm, n, p))
  p_score2 <- array(NA, dim = c(n_perm, n, n_pathway))
  y_pred <- rep(NA, n)
  valid_ind <- list()

  if (is.null(shuffle)) shuffle <- sample(n)
  if (is.null(active_var)) active_var <- seq(p)

  for (j in 1:k_fold){
    for(k in 1:2){
      if(k == 1){
        val.1 <- shuffle[floor(2*(j-1)*n/k_fold/2+1):floor((2*j-1)*n/k_fold/2)]
        val.2 <- shuffle[floor((2*j-1)*n/k_fold/2+1):floor(j*n/k_fold)]
      }else{
        val.2 <- shuffle[floor(2*(j-1)*n/k_fold/2+1):floor((2*j-1)*n/k_fold/2)]
        val.1 <- shuffle[floor((2*j-1)*n/k_fold/2+1):floor(j*n/k_fold)]
      }
      valid_ind[[2*(j-1)+k]] <- val.2
      hfs.object <- hybrid.corr(X[-c(val.1, val.2), ], y[-c(val.1, val.2)],
                                tau = 0.5, prop = prop, isAll = isAll,
                                utilities = utilities, M = M)
      feature.keep <- hfs.object$idx
      val_idx <- which(seq(dim(X)[1])[-val.2] %in% val.1)
      try(tau.object <- tune.tau(X=X[-val.2, ], y=y[-val.2],
            object=hfs.object,
            val_idx=val_idx, mod_fun=mod_fun,
            predict_fun=predict_fun, method=method,
            alpha=alpha, prop=prop, isAll=isAll, ...))
      try(feature.keep <- tau.object$idx)
      feature.keep <- feature.keep[feature.keep %in% active_var]

      if(hfs.object$type == "regression"){
        train_obj <- importDnnet(x=X[-val.2, feature.keep], y=y[-val.2])
        val_obj <- importDnnet(x=X[val.2, feature.keep], y=y[val.2])
      }else{
        train_obj <- importDnnet(x=X[-val.2, feature.keep], y=as.factor(y[-val.2]))
        val_obj <- importDnnet(x=X[val.2, feature.keep], y=as.factor(y[val.2]))
      }
      perm_mod <- permfit_inner(train=train_obj, validate = val_obj,
                                k_fold=0, n_perm=n_perm, mod_fun=mod_fun,
                                predict_fun = predict_fun,
                                pathway_list = pathway_list,
                                active_var = NULL, method=method,
                                shuffle = shuffle, ...)

      p_score2[,val.2,feature.keep] <- perm_mod$p_score2
      if (length(pathway_list) > 0) {
        p_score[,val.2,] <- perm_mod$p_score
      }
      y_pred[, val.2] <- perm_mod$y_pred
    }
  }
  
  if(is.null(colnames(X))) {
    imp <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    imp <- data.frame(var_name = colnames(X))
  }
  imp$importance <- apply(apply(p_score2, 2:3, mean), 2, mean, na.rm = TRUE)
  imp$importance_sd <- apply(apply(p_score2, 2:3, mean), 2,
                    function(x) sqrt(stats::var(x, na.rm=TRUE)/sum(!is.na(x))))
  imp$importance_pval <- 1 - stats::pnorm(imp$importance/imp$importance_sd)
  if(n_perm > 1) {
    imp$importance_sd_x <- apply(apply(p_score2, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    imp$importance_pval_x <- 1 - stats::pnorm(imp$importance/imp$importance_sd_x)
  }
  
  imp_block <- data.frame()
  if(n_pathway >= 1) {
    
    if(is.null(names(pathway_list))) {
      imp_block <- data.frame(block = paste0("P", 1:n_pathway))
    } else {
      imp_block <- data.frame(block = names(pathway_list))
    }
    imp_block$importance <- apply(apply(p_score, 2:3, mean), 2, mean, na.rm = TRUE)
    imp_blokc$importance_sd <- apply(apply(p_score2, 2:3, mean), 2,
                    function(x) sqrt(stats::var(x, na.rm=TRUE)/sum(!is.na(x))))
    imp_block$importance_pval <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd)
    if(n_perm > 1) {
      imp_block$importance_sd_x <- apply(apply(p_score, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
      imp_block$importance_pval_x <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd_x)
    }
  }
  return(new("PermFIT", model = NULL, importance = imp, block_importance = imp_block,
             validation_index = valid_ind, y_hat = y_pred))
}