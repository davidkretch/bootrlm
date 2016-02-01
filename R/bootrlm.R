#' Fit bootstrapped robust linear models
#'
#' @param formula A formula describing the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param r An integer number of replicates.
#' @param method A character string specifying the estimation method.
#' @param options A named list with options for particular estimation methods.
#' @return An object of class \code{bootrlm}, containing the following
#' components:
#' \item{r}{The integer number of replicates.}
#' \item{replicates}{A matrix whose columns are the row indices of each
#' replicate.}
#' \item{coefficients}{A matrix whose columns are the parameter
#' estimates for each replicate.}
#' \item{fitted}{A matrix whose columns are the fitted values
#' corresponding to the matrix \code{replicates}.}
#' \item{residuals}{A matrix whose columns are the residuals
#' corresponding to the matrix \code{replicates}.}
#' \item{scale}{A vector with the estimated scale for each replicate.}
#' \item{seed}{The seed used.}
#' \item{formula}{The formula supplied.}
#' \item{data}{The data frame supplied.}
#' \item{method}{The estimation method requested.}
#' @examples
#' library(MASS)
#' data(stackloss)
#' bootrlm_fit <- bootrlm(stack.loss ~ ., stackloss, r = 1000, method = "MM")
#' @export
bootrlm <- function(formula, data, r, method, options) {
  seed <- NA

  # construct data
  y <- model.response(model.frame(formula, data))
  x <- model.matrix(formula, data)
  n <- nrow(data)

  # construct replicates
  if (r > 0) {
    replicates <- matrix(sample(seq(0, n - 1), n * r, replace = TRUE), ncol = r)
  } else {
    replicates <- matrix(0:(n - 1), ncol = 1)
  }

  # determine options
  options <- list("k" = 1.548)

  # call estimation function
  fit <- bootrlm_cpp(y, x, replicates, 1, options)

  # format bootrlm object
  class(fit) <- "bootrlm"

  fit$replicates <- fit$replicates + 1 # indices start at 1
  fit$seed <- seed
  fit$formula <- formula
  fit$data <- data
  fit$method <- method

  return(fit)
}


#' @export
print.bootrlm <- function(x, ...) {
  str(x)
}
