###########################################################
### Define functions: importDnnet and splitDnnet

#' Import Data to create a \code{dnnetInput} object.
#'
#' @param x A \code{matrix} containing all samples/variables. It has to be \code{numeric}
#' and cannot be left blank. Any variable with missing value will be removed.
#' @param y A \code{numeric} or \code{factor} vector, indicating a continuous outcome or class label.
#' @param w A \code{numeric} vector, the sample weight. Will be 1 if left blank.
#'
#' @return An \code{dnnetInput} object.
#'
#' @importFrom methods new
#'
#' @seealso
#' \code{\link{dnnetInput-class}}
#' @export
importDnnet <- function(x, y, w = rep(1, dim(x)[1])) {

  new("dnnetInput", x = as.matrix(x), y = y, w = w)
}

#' A function to generate indice
#'
#' @param split As in \code{\link{dnnetInput-class}}.
#' @param n Sample size
#'
#' @return Returns a integer vector of indice.
#'
#' @seealso
#' \code{\link{dnnetInput-class}}
#'
#' @export
getSplitDnnet <- function(split, n) {

  if(is.numeric(split) && length(split) == 1 && split < 1)
    split <- sample(n, floor(n * split))

  if(is.numeric(split) && length(split) == 1 && split > 1)
    split <- 1:split

  if(is.character(split) && length(split) == 1 && split == "bootstrap")
    split <- sample(n, replace = TRUE)

  split
}

#' A function to split the \code{dnnetInput} object into a list of two \code{dnnetInput} objects:
#' one names train and the other named valid.
#'
#' @param object A \code{dnnetInput} object.
#' @param split A character, numeric variable or a numeric vector declaring a way to split
#' the \code{dnnetInput}. If it's number between 0 and 1, all samples will be split into two subsets
#' randomly, with the \code{train} containing such proportion of all samples and \code{valid} containing
#' the rest. If split is a character and is "bootstrap", the \code{train} will be a bootstrap sample
#' of the original data set and the \code{valid} will contain out-of-bag samples. If split is a vector
#' of integers, the \code{train} will contain samples whose indice are in the vector, and \code{valid} will
#' contain the rest.
#'
#' @return Returns a list of two \code{dnnetInput} objects.
#'
#' @seealso
#' \code{\link{dnnetInput-class}}
#'
#' @export
splitDnnet <-function(object, split) {

  split <- getSplitDnnet(split, dim(object@x)[1])

  train <- object
  train@x <- as.matrix(object@x[split, ])
  if(class(object@y)[1] != "matrix") {
    train@y <- object@y[split]
  } else {
    train@y <- object@y[split, ]
  }
  train@w <- object@w[split]
  if(class(object) == "dnnetSurvInput")
    train@e <- object@e[split]

  valid <- object
  valid@x <- as.matrix(object@x[-split, ])
  if(class(object@y)[1] != "matrix") {
    valid@y <- object@y[-split]
  } else {
    valid@y <- object@y[-split, ]
  }
  valid@w <- object@w[-split]
  if(class(object) == "dnnetSurvInput")
    valid@e <- object@e[-split]

  list(train = train, valid = valid, split = split)
}

###########################################################
### Model passed to PermFIT (Internal)

#' Model passed to PermFIT (Internal)
#'
#' Model passed to PermFIT (Internal)
#'
#' @param method Name of the model.
#' @param model.type Type of model.
#' @param object A dnnetInput object.
#' @param ... Orger parameters passed to the model.
#'
#' @return Returns a specific model.
#'
#' @importFrom randomForest randomForest
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom stats lm
#' @importFrom stats glm
#' @importFrom e1071 svm
#' @importFrom e1071 tune.svm
#' @importFrom xgboost xgboost
#'
#' @export
mod_permfit <- function(method, model.type, object, ...) {
  if (method == "random_forest") {

    mod <- do.call(randomForest::randomForest,
                     appendArg(appendArg(list(...), "x", object@x, TRUE), "y", object@y, TRUE))
  } else if (method == "lasso") {

    lasso_family <- ifelse(model.type == "regression", "gaussian","binomial")
    cv_lasso_mod <- glmnet::cv.glmnet(object@x, object@y, family = lasso_family)
    mod <- glmnet::glmnet(object@x, object@y, family = lasso_family,
                          lambda = cv_lasso_mod$lambda[which.min(cv_lasso_mod$cvm)])
  } else if (method == "linear") {

    if(model.type == "regression") {
      mod <- stats::lm(y ~ ., data.frame(x = object@x, y = object@y))
    } else if(model.type == "binary-classification") {
      mod <- stats::glm(y ~ ., family = "binomial", data = data.frame(x = object@x, y = object@y))
    }
  } else if (method == "svm") {

    if(model.type == "regression") {
      mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                             tunecontrol = e1071::tune.control(cross = 5))
      mod <- mod$best.model
    } else if(model.type == "binary-classification") {
      mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                             tunecontrol = e1071::tune.control(cross = 5))
      mod <- e1071::svm(object@x, object@y, gamma = mod$best.parameters$gamma, cost = mod$best.parameters$cost, probability = TRUE)
    } else {
      return("Not Applicable")
    }
  } else if (method == "xgboost") {

    if(model.type == "regression") {
      arg_xg <- list(...) %>%
        appendArg("data", object@x, TRUE) %>%
        appendArg("label", object@y, TRUE) %>%
        appendArg("verbose", 0, TRUE)
      if(!"nrounds" %in% names(arg_xg))
        arg_xg <- arg_xg %>% appendArg("nrounds", 50, TRUE)
      mod <- do.call(xgboost::xgboost, arg_xg)
    } else if(model.type == "binary-classification") {
      arg_xg <- list(...) %>%
        appendArg("data", object@x, TRUE) %>%
        appendArg("label", (object@y == levels(object@y)[1])*1, TRUE) %>%
        appendArg("verbose", 0, TRUE)
      if(!"nrounds" %in% names(arg_xg))
        arg_xg <- arg_xg %>% appendArg("nrounds", 50, TRUE)
      mod <- do.call(xgboost::xgboost, arg_xg)
    } else {
      return("Not Applicable")
    }
  } else {
    return("Not Applicable")
  }
  return(mod)
}

###########################################################
### Model prediction passed to PermFIT (Internal)

#' Model prediction passed to PermFIT (Internal)
#'
#' Model prediction passed to PermFIT (Internal)
#'
#' @param mod Model for prediction.
#' @param object A dnnetInput object.
#' @param method Name of the model.
#' @param model.type Type of the model.
#'
#' @return Returns predictions.
#'
#' @export
predict_mod_permfit <- function(mod, object, method, model.type) {

  if(model.type == "regression") {

    if(!method %in% c("linear", "lasso")) {
      return(predict(mod, object@x))
    } else if(method == "linear") {
      return(predict(mod, data.frame(x = object@x)))
    } else {
      return(predict(mod, object@x)[, "s0"])
    }
  } else if(model.type == "binary-classification") {

    if(method == "random_forest") {
      return(predict(mod, object@x, type = "prob")[, 1])
    } else if(method == "lasso") {
      return(1 - predict(mod, object@x, type = "response")[, "s0"])
    } else if (method == "linear") {
      return(1 - predict(mod, data.frame(x = object@x, y = object@y), type = "response"))
    } else if(method == "svm") {
      return(attr(predict(mod, object@x, decision.values = TRUE, probability = TRUE),
                  "probabilities")[, levels(object@y)[1]])
    } else if(method == "xgboost") {
      return(predict(mod, object@x))
    }
  } else {
    return("Not Applicable")
  }
}

###########################################################
### Log-likelihood Difference (Internal)

#' Log-likelihood Difference (Internal)
#'
#' Log-likelihood Difference (Internal)
#'
#' @param model.type Type of the model.
#' @param y_hat Y hat.
#' @param y_hat0 Another Y hat.
#' @param object Data object.
#' @param y_max Max prob for binary Y.
#' @param y_min Min prob for binary Y.
#'
#' @return Returns log-likelihood difference.
#'
#' @export
log_lik_diff <- function(model.type, y_hat, y_hat0, object, y_max = 1-10**-10, y_min = 10**-10) {

  if(model.type == "regression") {
    return((object@y - y_hat)**2 - (object@y - y_hat0)**2)
  } else if(model.type == "binary-classification") {
    y_hat <- ifelse(y_hat < y_min, y_min, ifelse(y_hat > y_max, y_max, y_hat))
    y_hat0 <- ifelse(y_hat0 < y_min, y_min, ifelse(y_hat0 > y_max, y_max, y_hat0))
    return(-(object@y == levels(object@y)[1])*log(y_hat) - (object@y != levels(object@y)[1])*log(1-y_hat) +
             (object@y == levels(object@y)[1])*log(y_hat0) + (object@y != levels(object@y)[1])*log(1-y_hat0))
  } else {
    return("Not Applicable")
  }
}

#' PermFIT: A permutation-based feature importance test. (Internal)
#'
#' @param train An dnnetInput object.
#' @param validate A validation dataset is required when k_fold = 0.
#' @param k_fold K-fold cross-fitting. If not, set k_fold to zero.
#' @param n_perm Number of permutations repeated.
#' @param mod_fun Function to train the model. Default \code{mod_permfit}.
#' @param predict_fun Function to predict the response. Default \code{predict_mod_permfit}.
#' @param pathway_list A list of pathways to be jointly tested.
#' @param active_var Index of active variables.
#' @param method Models, including \code{random_forest} for random forests,
#'  \code{lasso} for linear/logistic lasso, \code{linear} for
#'  linear/logistic regression,
#'   \code{svm} for svms with Gaussian kernels.
#' @param shuffle If shuffle is null, the data will be shuffled for
#'   cross-fitting; if random shuffle is not desired, please provide
#'   a bector of numbers for cross-fitting indices
#' @param ... Additional parameters passed to each method.
#'
#' @return Returns internal important score outputs.
#'
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats var
#'
#' @export
permfit_inner <- function(train, validate = NULL, k_fold = 5,
                    n_perm = 100, mod_fun=mod_permfit,
                    predict_fun=predict_mod_permfit,
                    pathway_list = list(),
                    active_var = NULL,
                    method = c("random_forest", "lasso",
                               "linear", "svm", "xgboost")[1],
                    shuffle = NULL,...) {
  n_pathway <- length(pathway_list)
  n <- dim(train@x)[1]
  if (is.null(active_var)){
    p <- dim(train@x)[2]
    active_var <- seq(p)
  }
  else{
    p <- length(active_var)
  }
  if(class(train) == "dnnetInput") {
    if(is.factor(train@y)) {
      model.type <- "binary-classification"
    } else {
      model.type <- "regression"
    }
  } else {
    stop("'train' has to be either a dnnetInput object.")
  }

  if(k_fold == 0) {

    if(is.null(validate))
      stop("A validation set is required when k = 0. ")
    n_valid <- dim(validate@x)[1]

    mod <- mod_fun(method, model.type, train, ...)
    f_hat_x <- predict_fun(mod, validate, method, model.type)
    valid_ind <- list(seq_along(validate@y))
    y_pred <- f_hat_x
    p_score <- array(NA, dim = c(n_perm, n_valid, 1))
    if(n_pathway >= 1) {
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      for(i in 1:n_pathway) {
        x_i <- validate@x
        for(l in 1:n_perm) {
          x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(n_valid), ]
          pred_i <- predict_fun(mod, importDnnet(x = x_i, y = validate@y), method, model.type)
          p_score[l, , i] <- log_lik_diff(model.type, pred_i, f_hat_x, validate)
        }
      }
    }

    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    j <- 1
    for(i in active_var) {
      x_i <- validate@x
      for(l in 1:n_perm) {
        x_i[, j] <- sample(x_i[, j])
        pred_i <- predict_fun(mod, importDnnet(x = x_i, y = validate@y), method, model.type)
        p_score2[l, , j] <- log_lik_diff(model.type, pred_i, f_hat_x, validate)
      }
      j <- j + 1
    }
  } else {
    valid_ind <- list()
    if(is.null(shuffle)) shuffle <- sample(n)
    n_valid <- n
    y_pred <- numeric(length(train@y))
    p_score <- array(NA, dim = c(n_perm, n_valid, 1))
    if(n_pathway >= 1)
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    valid_error <- numeric(k_fold)
    for(k in 1:k_fold) {
      train_spl <- splitDnnet(train, shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)])
      valid_ind[[k]] <- shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)]

      mod <- mod_fun(method, model.type, train_spl$valid, ...)
      f_hat_x <- predict_fun(mod, train_spl$train, method, model.type)
      valid_error[k] <- sum(log_lik_diff(model.type, f_hat_x, f_hat_x, train_spl$train))
      y_pred[valid_ind[[k]]] <- f_hat_x
      if(k == 1) {

        final_model <- mod
      }

      if(n_pathway >= 1) {
        for(i in 1:n_pathway) {
          for(l in 1:n_perm) {
            x_i <- train_spl$train@x
            x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(dim(x_i)[1]), ]
            pred_i <- predict_fun(mod, importDnnet(x = x_i, y = train_spl$train@y), method, model.type)
            p_score[l, shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] <- log_lik_diff(model.type, pred_i, f_hat_x, train_spl$train)
          }
        }
      }
      j <- 1
      for(i in active_var) {
        x_i <- train_spl$train@x
        for(l in 1:n_perm) {
          x_i[, j] <- sample(x_i[, j])
          pred_i <- predict_fun(mod, importDnnet(x = x_i, y = train_spl$train@y), method, model.type)
          p_score2[l, shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], j] <- log_lik_diff(model.type, pred_i, f_hat_x, train_spl$train)
        }
        j <- j + 1
      }
    }
    mod <- final_model
    valid_error <- sum(valid_error)/n_valid
  }

  return(list(model = mod, p_score = p_score, p_score2 = p_score2, n_pathway=n_pathway,
            valid_ind=valid_ind, y_pred=y_pred, n_valid=n_valid))
}

###########################################################
### PermFIT

#' PermFIT: A permutation-based feature importance test.
#'
#' @param train An dnnetInput object.
#' @param validate A validation dataset is required when k_fold = 0.
#' @param k_fold K-fold cross-fitting. If not, set k_fold to zero.
#' @param n_perm Number of permutations repeated.
#' @param mod_fun Function to train the model. Default \code{mod_permfit}.
#' @param predict_fun Function to predict the response. Default \code{predict_mod_permfit}.
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
#' @param ... Additional parameters passed to each method.
#'
#' @return Returns a PermFIT object.
#'
#' @importFrom stats sd
#' @importFrom stats pnorm
#' @importFrom stats var
#'
#' @export
permfit <- function(train, validate = NULL, k_fold = 5,
                    n_perm = 100, mod_fun = mod_permfit,
                    predict_fun = predict_mod_permfit,
                    pathway_list = list(),
                    active_var = NULL,
                    method = c("random_forest", "lasso",
                               "linear", "svm", "xgboost")[1],
                    shuffle = NULL, ...) {
  if (is.null(active_var)){
    active_var <- seq(dim(train@x)[2])
  }

  inner <- permfit_inner(train, validate, k_fold, n_perm, mod_fun,
                         predict_fun, pathway_list, active_var,
                         method, shuffle, ...)

  mod <- inner$mod
  p_score <- inner$p_score
  p_score2 <- inner$p_score2
  n_pathway <- inner$n_pathway
  valid_ind <- inner$valid_ind
  y_pred <- inner$y_pred
  n_valid <- inner$n_valid
  remove(inner)

  if(is.null(colnames(train@x))) {
    imp <- data.frame(var_name = paste0("V", active_var))
  } else  {
    imp <- data.frame(var_name = colnames(train@x)[active_var])
  }
  imp$importance <- apply(apply(p_score2, 2:3, mean), 2, mean, na.rm = TRUE)
  imp$importance_sd <- sqrt(apply(apply(p_score2, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
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
    imp_block$importance_sd <- sqrt(apply(apply(p_score, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
    imp_block$importance_pval <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd)
    if(n_perm > 1) {
      imp_block$importance_sd_x <- apply(apply(p_score, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
      imp_block$importance_pval_x <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd_x)
    }
  }

  return(new("PermFIT", model = mod, importance = imp, block_importance = imp_block,
             validation_index = valid_ind, y_hat = y_pred))
}
