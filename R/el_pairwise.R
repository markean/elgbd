#' Pairwise comparisons for general block designs with empirical likelihood
#'
#' @description
#'   Tests all pairwise comparisons or comparisons with control for general
#'   block designs with empirical likelihood. Two single step asymptotic
#'   \eqn{k}-FWER (generalized family-wise error rate) controlling procedures
#'   are available: asymptotic Monte Carlo (AMC) and nonparametric bootstrap
#'   (NB).
#'
#' @param formula
#'   An object of class [`formula`] (or one that can be coerced to that class)
#'   for a symbolic description of the model to be fitted. It must specify the
#'   variables for response, treatment, and block as `response ~ treatment |
#'   block`. Note that the use of vertical bar (`|`) separating treatment and
#'   block.
#' @param data
#'   A data frame, list or environment (or object coercible by [as.data.frame()]
#'   to a data frame) containing the variables in `formula`.
#' @param control
#'   An optional single character that specifies the treatment for comparisons
#'   with control.
#' @param k
#'   A single integer for \eqn{k} in \eqn{k}-FWER. Defaults to `1`.
#' @param alpha
#'   A single numeric for the overall significance level. Defaults to `0.05`.
#' @param method
#'   A single character for the procedure to be used; either `AMC` or `NB` is
#'   supported. Defaults to `AMC`.
#' @param B
#'   A single integer for the number of Monte Carlo samples for the AMC (number
#'   of bootstrap replicates for the NB).
#' @param nthreads
#'   A single integer for the number of threads for parallel computation via
#'   'OpenMP' (if available). Defaults to `1`.
#' @param maxit
#'   A single integer for the maximum number of iterations for constrained
#'   minimization of empirical likelihood. Defaults to `10000`.
#' @param abstol
#'   A single numeric for the the absolute convergence tolerance for
#'   optimization. Defaults to `1e-08`.
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
#' # All pairwise comparisons
#' data("clothianidin")
#' el_pairwise(clo ~ trt | blk, data = clothianidin, B = 1000)
#'
#' # Comparisons with control
#' el_pairwise(clo ~ trt | blk,
#'   control = "Naked", data = clothianidin, B = 1000
#' )
#' @export
el_pairwise <- function(formula, data, control = NULL, k = 1L, alpha = 0.05,
                        method = c("AMC", "NB"), B, nthreads = 1L,
                        maxit = 10000L, abstol = 1e-08, verbose = FALSE) {
  alpha <- validate_alpha(alpha)
  B <- validate_b(B)
  max_threads <- get_max_threads()
  nthreads <- validate_nthreads(nthreads, max_threads)
  maxit <- validate_maxit(maxit)
  abstol <- validate_tol(abstol)
  method <- match.arg(method)
  f <- attributes(terms(formula))
  stopifnot(
    "`formula` must be specified as `response ~ treatment | block`." =
      (isTRUE(f$response == 1L && f$intercept == 1L)),
    "`formula` must be specified as `response ~ treatment | block`." =
      (isTRUE(length(f$variables) == 3L)),
    "`formula` must be specified as `response ~ treatment | block`." =
      (isTRUE(typeof(f$variables[[3L]]) == "language" &&
        length(f$variables[[3L]]) == 3L)),
    "`|` operator is missing." = (isTRUE(f$variables[[3L]][[1L]] == "|")),
    "Transformation on variables is not allowed." =
      (isTRUE(typeof(f$variables[[3L]][[2L]]) == "symbol" &&
        typeof(f$variables[[3L]][[3L]]) == "symbol")),
    "Specify distinct variables for the treatments and blocks." =
      (isTRUE(f$variables[[3L]][[2L]] != f$variables[[3L]][[3L]]))
  )
  # Pseudo formula for `model.frame`
  lhs <- f$variables[[2L]]
  rhs <- c(f$variables[[3L]][[2L]], f$variables[[3L]][[3L]])
  pf <- formula(paste(lhs, paste(rhs, collapse = " + "), sep = " ~ "))
  # Extract `model.frame`
  mf <- match.call()
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf[[2L]] <- pf
  mf <- eval(mf, parent.frame())
  attributes(mf)$terms <- NULL
  ## Type conversion
  mf[[1L]] <- as.numeric(mf[[1L]])
  mf[[2L]] <- as.factor(mf[[2L]])
  mf[[3L]] <- as.factor(mf[[3L]])
  if (nlevels(mf[[2L]]) >= nlevels(mf[[3L]])) {
    stop("The number of blocks must be larger than the number of treatments.")
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
  # Replace `NA` with `0`
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
  # Check whether all pairwise comparisons or comparisons to control
  match.arg(control, gbd$trt)
  if (is.null(control)) {
    ctrl <- 0
  } else {
    ctrl <- which(control == gbd$trt)
  }

  # Pairwise comparisons
  out <- pairwise(gbd$model_matrix, gbd$incidence_matrix,
    control = ctrl, k, alpha, interval = TRUE,
    method, B, nthreads, threshold = 50, maxit, abstol
  )
  out$trt <- gbd$trt
  out$control <- control
  out$model.matrix <- gbd$model_matrix
  out$incidence.matrix <- gbd$incidence_matrix
  if (isTRUE(verbose)) {
    if (!all(out$convergence)) {
      if (method == "NB") {
        message(
          "Convergence failed and switched to AMC for confidence intervals.\n"
        )
      } else {
        message("Convergence failed.\n")
      }
    }
  }
  out
}
