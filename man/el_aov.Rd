% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/el_aov.R
\name{el_aov}
\alias{el_aov}
\title{Fit an Analysis of Variance Model with Empirical Likelihood}
\usage{
el_aov(formula, data, maxit = 10000, abstol = 1e-08)
}
\arguments{
\item{formula}{A formula object. It must specify variables for response and treatment as 'response ~ treatment'.}

\item{data}{A data frame containing the variables in the formula.}

\item{maxit}{Maximum number of iterations for optimization. Defaults to 10000.}

\item{abstol}{Absolute convergence tolerance for optimization. Defaults to 1e-08.}
}
\value{
A list with class \code{c("el_aov", "melt")}.
}
\description{
Fit an one-way analysis of variance model with empirical likelihood.
}
\examples{
data("clothianidin")
el_aov(clo ~ trt, clothianidin)
}
\references{
Owen, Art. 1991. “Empirical Likelihood for Linear Models.” The Annals of Statistics 19 (4). \doi{10.1214/aos/1176348368}.
}
\seealso{
\link{el_test}
}
