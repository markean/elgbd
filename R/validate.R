#' Validate `alpha`
#'
#' Validate `alpha` in [elt()].
#'
#' @param alpha A single numeric.
#' @return A single numeric.
#' @noRd
validate_alpha <- function(alpha) {
  stopifnot(
    "`alpha` must be a finite single numeric." =
      isTRUE(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha)),
    "`alpha` must be between 0 and 1." = isTRUE(alpha > 0 && alpha < 1)
  )
  alpha
}

#' Validate `b`
#'
#' Validate `b` in [el_control()].
#'
#' @param b A single integer.
#' @return A single integer.
#' @noRd
validate_b <- function(b) {
  b <- tryCatch(as.integer(b), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`b` must be a finite single integer." = isTRUE(!is.na(b)),
    "`b` must be a positive single integer." = b > 0L
  )
  b
}

#' Validate `nthreads`
#'
#' Validate `nthreads` in [el_control()].
#'
#' @param nthreads A single integer.
#' @param max_threads A single integer.
#' @return A single integer.
#' @noRd
validate_nthreads <- function(nthreads, max_threads) {
  nthreads <- tryCatch(as.integer(nthreads), warning = \(x) NA, error = \(x) NA)
  stopifnot("`nthreads` must be a single integer." = isTRUE(!is.na(nthreads)))
  if (nthreads < 1) {
    warning("`nthreads` is set to 1.")
    nthreads <- 1L
  }
  if (nthreads > max_threads) {
    warning("`nthreads` is set to the maximum number of threads available.")
    nthreads <- max_threads
  }
  nthreads
}

#' Validate `maxit`
#'
#' Validate `maxit` in [el_control()].
#'
#' @param maxit A single integer.
#' @return A single integer.
#' @noRd
validate_maxit <- function(maxit) {
  maxit <- tryCatch(as.integer(maxit), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`maxit` must be a finite single integer." = isTRUE(!is.na(maxit)),
    "`maxit` must be a positive single integer." = maxit > 0L
  )
  maxit
}

#' Validate `tol`
#'
#' Validate `tol` in [el_control()].
#'
#' @param tol A single numeric.
#' @return A single numeric.
#' @noRd
validate_tol <- function(tol) {
  tol <- tryCatch(as.numeric(tol), warning = \(x) NA, error = \(x) NA)
  stopifnot(
    "`tol` must be a finite single numeric." =
      isTRUE(!is.na(tol) && is.finite(tol)),
    "`tol` is too small." = tol >= .Machine$double.eps
  )
  tol
}
