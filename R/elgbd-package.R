#' @references
#'   Kim E, MacEachern SN, Peruggia M (2023).
#'   "Empirical likelihood for the analysis of experimental designs."
#'   \emph{Journal of Nonparametric Statistics}, **35**(4), 709--732.
#'   \doi{10.1080/10485252.2023.2206919}.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom stats .getXlevels
#' @importFrom stats coef
#' @importFrom stats complete.cases
#' @importFrom stats printCoefmat
#' @importFrom stats qchisq
#' @importFrom stats reshape
#' @importFrom stats setNames
#' @importFrom stats terms
#' @importFrom stats weights
#' @useDynLib elgbd, .registration = TRUE
## usethis namespace: end
NULL
