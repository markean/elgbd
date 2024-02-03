#include "utils_gbd.h"

// [[Rcpp::export]]
Rcpp::List ELtest(const Eigen::MatrixXd& x,
                  const Eigen::MatrixXd& c,
                  const Eigen::MatrixXd& lhs,
                  const Eigen::VectorXd& rhs,
                  const double threshold,
                  const int maxit = 1e4,
                  const double abstol = 1e-8) {
  // Check lhs & rhs
  if (lhs.rows() > lhs.cols()) {
    Rcpp::stop("nrow(lhs) must not exceed ncol(lhs)");
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(lhs);
  if (lu_decomp.rank() != lhs.rows()) {
    Rcpp::stop("lhs must have full rank");
  }
  if (lhs.rows() != rhs.rows()) {
    Rcpp::stop("dimensions of lhs and rhs do not match");
  }

  // Global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();

  minEL el =
    test_gbd_EL(theta_hat, x, c, lhs, rhs, threshold, maxit, abstol);

  Rcpp::List result;
  result["coefficients"] = theta_hat;
  result["statistic"] = 2 * el.nlogLR;
  result["optim"] = Rcpp::List::create(
    Rcpp::Named("par") = el.par,
    Rcpp::Named("lambda") = el.lambda,
    Rcpp::Named("iterations") = el.iterations,
    Rcpp::Named("convergence") = el.convergence);
  result.attr("class") = "elgbd";
  return result;
}
