#' @export
print.el_aov <- function(x, ...) {
  stopifnot(inherits(x, "elgbd"))
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nMinimizer:\n")
  cat(format(round(x$optim$par, 4L), scientific = FALSE))
  cat("\n\n")
  cat("Statistic:\n")
  cat(format(round(x$statistic, 4L), scientific = FALSE))
  cat("\n\n")
}

#' @export
print.pairwise <- function(x, ...) {
  stopifnot(inherits(x, "elgbd"))
  cat("\n\tEmpirical Likelihood Multiple Tests\n\n")
  if (is.null(x$control)) {
    cat("All pairwise comparisons\n\n")
    rname <- vector("character", length = 0L)
    for (i in 1L:(length(x$trt) - 1L)) {
      for (j in (i + 1L):length(x$trt)) {
        rname <- c(rname, paste(x$trt[i], "-", x$trt[j]))
      }
    }
  } else {
    cat("Comparisons with control\n\n")
    diff <- setdiff(x$trt, x$control)
    rname <- vector("character", length = length(diff))
    for (i in seq_along(diff)) {
      rname[i] <- paste(diff[i], "-", x$control)
    }
  }
  out <- data.frame(
    row.names = rname, Estimate = x$estimate, Chisq = x$statistic,
    Lwr.ci = x$lower, Upr.ci = x$upper, p.adj = x$p.adj
  )
  printCoefmat(out,
    digits = min(4L, getOption("digits")),
    dig.tst = min(3L, getOption("digits")), cs.ind = c(1L, 3L, 4L),
    tst.ind = 2L, P.values = TRUE, has.Pvalue = TRUE, eps.Pvalue = 1e-03
  )
  cat("\n")
  cat(paste(c("k", "level", "method", "cutoff"),
    c(x$k, x$level, x$method, round(x$cutoff, 4L)),
    collapse = ", ", sep = ": "
  ))
  cat("\n\n")
}

#' @export
print.el_test <- function(x, digits = getOption("digits"), ...) {
  stopifnot(inherits(x, "elgbd"))
  cat("\n\tEmpirical Likelihood Test\n\n")
  cat("General block designs\n")
  if (!is.null(x$coefficients)) {
    cat("\nMaximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  out <- character()
  if (!is.null(x$statistic)) {
    cat("\n")
    out <- c(out, paste(
      "Statistic: ", names(x$statistic),
      format(x$statistic, digits = max(1L, digits - 2L))
    ))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("\n")
  invisible(x)
}
