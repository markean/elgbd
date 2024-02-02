#' Analysis of variance with empirical likelihood
#'
#' @description
#'   Fits an one-way analysis of variance model with empirical likelihood.
#'
#' @param formula
#'   An object of class [`formula`] (or one that can be coerced to that class)
#'   for a symbolic description of the model to be fitted. It must specify
#'   the variables for response and treatment as `response ~ treatment`.
#' @param data
#'   A data frame containing the variables in `formula`.
#' @param maxit
#'   A single integer for the maximum number of iterations for optimization.
#'   Defaults to `10000`.
#' @param abstol
#'   A single numeric for the absolute convergence tolerance for optimization.
#'   Defaults to `1e-08`.
#' @return
#'   A list containing the model fit and optimization results.
#' @references
#'   Owen, A (1991).
#'   "Empirical Likelihood for Linear Models."
#'   \emph{The Annals of Statistics}, **19**(4), 1725--1747.
#'   \doi{10.1214/aos/1176348368}.
#' @examples
#' data("clothianidin")
#' el_aov(clo ~ trt, clothianidin)
#' @export
el_aov <- function(formula, data, maxit = 1e+04, abstol = 1e-08) {
  # Check formula
  f <- attributes(terms(formula))
  if (any(
    # Response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3,
    # No other formula
    typeof(f$variables[[3]]) != "symbol" ||
      length(f$variables[[3]]) != 1
  )
  ) {
    stop("invalied model formula. specify formula as 'response ~ treatment'")
  }

  # Extract model frame
  mf <- cl <- match.call()
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())

  # Type conversion
  # Response
  mf[[1L]] <- as.numeric(mf[[1L]])
  # Treatment
  mf[[2L]] <- as.factor(mf[[2L]])

  # Extract model terms & levels
  # Model terms
  mt <- attr(mf, "terms")
  # Number of levels
  p <- nlevels(mf[[2L]])
  if (p < 2L) {
    stop("contrasts can be applied only to factors with 2 or more levels")
  }
  # Levels
  lv <- .getXlevels(mt, mf)
  # Name for coefficients
  nm <- paste0(names(lv), lv[[1L]])

  # Construct a general block design
  # Incidence matrix
  c <- unclass(table(
    factor(row.names(mf), levels = unique(row.names(mf))),
    mf[[2L]]
  ))
  # Model matrix
  x <- mf[[1L]] * c

  # Specify hypothesis
  lhs <- matrix(0, nrow = p - 1, ncol = p)
  if (p == 2L) {
    lhs <- matrix(c(1, -1), nrow = 1L)
  } else {
    diag(lhs) <- 1
    diag(lhs[, -1L]) <- -1
  }
  rhs <- rep(0, p - 1)

  # Test hypothesis
  out <- ELtest(x, c, lhs, rhs, threshold = 500, maxit, abstol)
  out$coefficients <- setNames(out$coefficients, nm)
  out$xlevels <- lv
  out$call <- cl
  out$terms <- mt
  out$model <- list(model.matrix = x, incidence.matrix = c)
  class(out) <- c("el_aov", oldClass(out))
  if (!out$optim$convergence) {
    warning("convergence failed\n")
  }
  out
}
