#' Hypothesis testing with empirical likelihood
#'
#' @description
#'   Tests single hypothesis for general block designs with empirical
#'   likelihood.
#'
#' @param formula
#'   An object of class [`formula`] (or one that can be coerced to that class)
#'   for a symbolic description of the model to be fitted. It must specify the
#'   variables for response, treatment, and block as `response ~ treatment |
#'   block`. Note that the use of vertical bar (`|`) separating treatment and
#'   block.
#' @param data
#'   A data frame containing the variables in `formula`.
#' @param lhs
#'   A numeric matrix specifying the left-hand side of a hypothesis in terms of
#'   parameters.
#' @param rhs
#'   An optional numeric vector specifying the right-hand side the hypothesis.
#'   If not specified, it is set to the zero vector. Defaults to `NULL`.
#' @param maxit
#'   A single integer for the maximum number of iterations for optimization.
#'   Defaults to `10000`.
#' @param abstol
#'   A single numeric for the absolute convergence tolerance for optimization.
#'   Defaults to `1e-08`.
#' @param verbose
#'   A single logical. If `TRUE`, a message on the convergence status is
#'   printed. Defaults to `FALSE`.
#' @return
#'   A list containing the model fit and optimization results.
#' @references
#'   Kim E, MacEachern SN, Peruggia M (2023).
#'   "Empirical likelihood for the analysis of experimental designs."
#'   \emph{Journal of Nonparametric Statistics}, **35**(4), 709--732.
#'   \doi{10.1080/10485252.2023.2206919}.
#' @examples
#' # Test for equal means
#' data("clothianidin")
#' el_test(clo ~ trt | blk, clothianidin,
#'   lhs = matrix(c(
#'     1, -1, 0, 0,
#'     0, 1, -1, 0,
#'     0, 0, 1, -1
#'   ), byrow = TRUE, nrow = 3L)
#' )
#' @export
el_test <- function(formula, data, lhs, rhs = NULL, maxit = 1e+04,
                    abstol = 1e-08, verbose = FALSE) {
  # Check formula
  f <- attributes(terms(formula))
  if (any(
    # Response required & no arbitrary manipulation on intercept
    f$response == 0, f$intercept == 0,
    length(f$variables) != 3L,
    # No other formula
    typeof(f$variables[[3L]]) != "language" ||
      length(f$variables[[3L]]) != 3L,
    # "|" operator needed
    f$variables[[3L]][[1L]] != "|",
    # No transformation of variables
    typeof(f$variables[[3L]][[2L]]) != "symbol" ||
      typeof(f$variables[[3L]][[3L]]) != "symbol",
    # Distinct variables for treatment and block
    f$variables[[3L]][[2L]] == f$variables[[3L]][[3L]]
  )
  ) {
    stop("Invalied model formula. specify formula as `response ~ treatment | block`")
  }

  ## Pseudo formula for model frame
  l <- f$variables[[2L]]
  r <- c(f$variables[[3L]][[2L]], f$variables[[3L]][[3L]])
  pf <- formula(paste(l, paste(r, collapse = " + "), sep = " ~ "))

  ## Extract model frame
  mf <- match.call()
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf[[2L]] <- pf
  mf <- eval(mf, parent.frame())
  attributes(mf)$terms <- NULL

  # Type conversion
  # Response
  mf[[1L]] <- as.numeric(mf[[1L]])
  # Treatment
  mf[[2L]] <- as.factor(mf[[2L]])
  # block
  mf[[3L]] <- as.factor(mf[[3L]])
  if (nlevels(mf[[2L]]) >= nlevels(mf[[3L]])) {
    stop("number of blocks should be larger than number of treatments")
  }

  # Construct general block design
  # Incidence matrix
  c <- unclass(table(mf[[3L]], mf[[2L]]))
  # Model matrix
  x <- reshape(mf[order(mf[[2L]]), ],
    idvar = names(mf)[3L],
    timevar = names(mf)[2L],
    v.names = names(mf)[1L],
    direction = "wide"
  )
  x <- x[order(x[[names(mf)[3L]]]), ]
  # Replace NA with 0
  x[is.na(x)] <- 0
  # Remove block variable and convert to matrix
  x[names(mf)[3L]] <- NULL
  x <- as.matrix(x)
  # Name rows and columns
  dimnames(x) <- list(levels(mf[[3L]]), levels(mf[[2L]]))
  # General block design
  gbd <-
    list("model_matrix" = x, "incidence_matrix" = c, "trt" = levels(mf[[2L]]))
  class(gbd) <- c("gbd", "elgbd")

  # Test for lhs and rhs
  if (is.null(rhs)) {
    rhs <- rep(0, nrow(lhs))
  }

  # Test hypothesis
  out <- ELtest(gbd$model_matrix, gbd$incidence_matrix, lhs, rhs,
    threshold = nrow(lhs) * 500, maxit, abstol
  )
  out$trt <- gbd$trt
  out$model.matrix <- gbd$model_matrix
  out$incidence.matrix <- gbd$incidence_matrix
  class(out) <- c("el_test", oldClass(out))
  if (isTRUE(verbose)) {
    if (!out$optim$convergence) {
      message("Convergence failed.\n")
    }
  }
  out
}
