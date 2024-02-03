#' @export
print.el_aov <- function(x, ...) {
  stopifnot(inherits(x, "elgbd"))
  cat("Call:\n")
  dput(x$call, control = NULL)
  cat("\nMinimizer:\n")
  cat(format(round(x$optim$par, 4L), scientific = FALSE))
  cat("\n\n")
  cat("Statistic:\n")
  cat(format(round(x$optim$n2logLR, 4L), scientific = FALSE))
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
  cat("\n")
  cat("Empirical Likelihood Test:", x$optim$type, "\n")
  cat("\n")
  out <- character()
  if (!is.null(x$statistic)) {
    out <- c(out, paste(
      "Chisq", names(x$statistic), "=",
      format(x$statistic, digits = max(1L, digits - 2L))
    ))
  }
  if (!is.null(x$df)) {
    out <- c(out, paste("df", "=", x$df))
  }
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value", if (startsWith(fp, "<")) {
      fp
    } else {
      paste(
        "=",
        fp
      )
    }))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
      if (length(x$null.value) == 1L) {
        alt.char <- switch(x$alternative,
          two.sided = "not equal to",
          less = "less than",
          greater = "greater than"
        )
        cat("true ", names(x$null.value), " is ", alt.char,
          " ", x$null.value, "\n",
          sep = ""
        )
      } else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits = digits, ...)
      }
    } else {
      cat(x$alternative, "\n", sep = "")
    }
  }
  if (!is.null(x$coefficients)) {
    cat("maximum EL estimates:\n")
    print(x$coefficients, digits = digits, ...)
  }
  cat("\n")
  invisible(x)
}
