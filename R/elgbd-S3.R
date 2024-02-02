#' @export
print.el_aov <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nminimizer:\n")
  cat(format(round(x$optim$par, 4), scientific = FALSE))
  cat("\n\n")
  cat("statistic:\n")
  cat(format(round(x$optim$n2logLR, 4), scientific = FALSE))
  cat("\n\n")
}

#' @export
print.pairwise <- function(x, ...) {
  stopifnot(inherits(x, "melt"))
  cat("\n")
  cat("Empirical Likelihood Multiple Hypothesis Testing\n\n")
  # set row names
  if (is.null(x$control)) {
    cat("Test: all pairwise comparisons\n\n")
    rname <- vector("character", length = 0)
    for (i in 1L:(length(x$trt) - 1L)) {
      for (j in (i + 1L):length(x$trt)) {
        rname <- c(rname, paste(x$trt[i], "-", x$trt[j]))
      }
    }
  } else {
    cat("Test: comparisons with control\n\n")
    diff <- setdiff(x$trt, x$control)
    rname <- vector("character", length = length(diff))
    for (i in 1L:length(diff)) {
      rname[i] <- paste(diff[i], "-", x$control)
    }
  }
  out <- data.frame(
    row.names = rname, estimate = x$estimate,
    statistic = x$statistic, lwr.ci = x$lower,
    upr.ci = x$upper,
    p.adj = round(x$p.adj, 4)
  )
  printCoefmat(out,
    digits = min(4L, getOption("digits")), cs.ind = c(1, 3, 4),
    tst.ind = 2L, dig.tst = min(3L, getOption("digits")),
    P.values = TRUE, has.Pvalue = TRUE, eps.Pvalue = 1e-03
  )
  cat("\n")
  cat(paste(c("k", "level", "method", "cutoff"),
    c(x$k, x$level, x$method, round(x$cutoff, 4)),
    collapse = ", ", sep = ": "
  ))
  cat("\n\n")
}
