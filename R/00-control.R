#################################
#' Append an argument on a list
#'
#' Append an argument on a list of arguments
#'
#' @param args The original argument list.
#' @param argName The name of the argument to add.
#' @param argValue The value of the argument to add.
#' @param forced a \code{logical} value indicating if the argument with the same
#'  name already existed, whether it should be added forcedly.
#'
#' @return An argument list.
#'
#' @export
appendArg <- function(args, argName, argValue, forced) {
  
  if((!argName %in% names(args)) && (!forced)) {
    
    # cat("Default behavior: setting", argName, "to", argValue, "...\n")
    args[[argName]] <- argValue
  } else if(forced) {
    
    # cat("Forced behavior: setting", argName, "to", argValue, "(cannot override) ...\n")
    args[[argName]] <- argValue
  }
  
  return(args)
}

#################################
#' Remove an argument on a list
#'
#' Remove an argument on a list of arguments
#'
#' @param args The original argument list.
#' @param argName The name of the argument to be remmoved. If the name does not exist,
#'  no argument will be removed.
#'
#' @return An argument list.
#'
#' @export
removeArg <- function(args, argName) {
  
  if(argName %in% names(args))
    args[[which(names(args) == argName)]] <- NULL
  return(args)
}

#################################
#### dnnetInput class
#' An S4 class containing predictors (x), response (y) and sample weights (w)
#'
#' @slot x A numeric matrix, the predictors
#' @slot y A factor or numeric vector, either the class labels or continuous responses
#' @slot w A numeric vector, sample weights
#'
#' @seealso
#' \code{\link{dnnet-class}}\cr
#' @export
setClass("dnnetInput",
         slots = list(
           x = "matrix",
           y = "ANY",
           w = "numeric"
         ))

#################################
#### PermFIT class
#' An S4 class containing an permutation-based feature importance test (PermFIT) object
#'
#' @slot model a model for ITR prediction
#' @slot importance a data.frame of importance
#' @slot block_importance a data.frame of block importance
#' @slot validation_index indices for cross-fitting
#' @slot y_hat cross-fitted prediction or prediction from validation
#'
#' @export
setClass("PermFIT",
         slots = list(
           model = "ANY",
           importance = "data.frame",
           block_importance = "data.frame",
           validation_index = "ANY",
           y_hat = "ANY"
         ))

#################################
#### dnnetInput class

#' @name show
#' @rdname show
#'
#' @title Method show for the package
#'
#' @description The method show for \code{dnnetInput} or \code{dnnet} object.
#'
#' @param object A \code{dnnet} or \code{dnnetInput} object.
#'
#' @seealso
#' \code{\link{dnnetInput-class}}\cr
#' \code{\link{dnnet-class}}\cr
#'
NULL

#' @rdname show
#' @importFrom methods show
#'
#' @export
setMethod("show",
          "dnnetInput",
          function(object) {

            print("y: (first 6)")
            if(class(object@y) == "matrix")
              print(object@y[1:min(6, dim(object@x)[1]), ])
            else
              print(object@y[1:min(6, dim(object@x)[1])])
            print("x: (first 6x6)")
            print(paste("dim(x) =", dim(object@x)[1], "x", dim(object@x)[2]))
            print(object@x[1:min(6, dim(object@x)[1]), 1:min(6, dim(object@x)[2])])
            print("w: (first 6)")
            print(object@w[1:min(6, dim(object@x)[1])])
          })