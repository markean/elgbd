#include "utils_pairwise.h"

// [[Rcpp::export]]
Rcpp::List pairwise(const Eigen::MatrixXd& x,
                    const Eigen::MatrixXd& c,
                    const int control = 0,
                    const int k = 1,
                    const double level = 0.05,
                    const bool interval = true,
                    const std::string method = "AMC",
                    const int B = 1e4,
                    const int nthread = 1,
                    const bool progress = true,
                    const double threshold = 50,
                    const int maxit = 1e4,
                    const double abstol = 1e-8) {
  if (level <= 0 || level >= 1) {
    Rcpp::stop("`level` must be between 0 and 1.");
  }
  // Pairs
  std::vector<std::array<int, 2>> pairs = comparison_pairs(x.cols(), control);
  // Number of hypotheses
  const int m = pairs.size();
  // Estimates
  std::vector<double> estimate(m);
  // Statistics
  std::vector<double> statistic(m);
  // Convergences
  std::vector<bool> convergence(m);

  // Global minimizer
  const Eigen::VectorXd theta_hat =
    x.array().colwise().sum() / c.array().colwise().sum();
  if (progress) {
    REprintf("Computing statistics...");
  }
  // 1. statistics
  for (int i = 0; i < m; ++i) {
    Rcpp::checkUserInterrupt();
    estimate[i] = theta_hat(pairs[i][0]) - theta_hat(pairs[i][1]);
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
    lhs(pairs[i][0]) = 1;
    lhs(pairs[i][1]) = -1;
    minEL pairwise_result =
      test_gbd_EL(theta_hat, x, c, lhs, Eigen::Matrix<double, 1, 1>(0),
                  threshold, maxit, abstol);
    statistic[i] = 2 * pairwise_result.nlogLR;
    convergence[i] = pairwise_result.convergence;
  }
  // If any of the statistics is not converged, switch
  bool anyfail = std::any_of(convergence.begin(), convergence.end(),
                             [](bool v) {return !v;});
  // Bootstrap statistics
  if (progress) {
    REprintf("\nComputing cutoff...");
  }
  Eigen::ArrayXd bootstrap_statistics_pairwise(B);
  if (method == "AMC") {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
  } else {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_NB(x, c, k, pairs, B,
                                       level, nthread, progress,
                                       threshold, maxit, abstol);
  }
  // 2. Adjusted p-values
  std::vector<double> adj_pvalues(m);
  for (int i = 0; i < m; ++i) {
    adj_pvalues[i] =
      static_cast<double>(
        (bootstrap_statistics_pairwise >= statistic[i]).count()) / B;
  }
  // 3. Cutoff
  Rcpp::Function quantile("quantile");
  double cutoff =
    Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                              Rcpp::Named("probs") = 1 - level));

  // Result
  Rcpp::List result;
  result["estimate"] = estimate;
  result["statistic"] = statistic;
  result["convergence"] = convergence;
  result["cutoff"] = cutoff;
  if (method == "NB" && anyfail) {
    bootstrap_statistics_pairwise =
      bootstrap_statistics_pairwise_AMC(x, c, k, pairs, B, level);
    cutoff =
      Rcpp::as<double>(quantile(bootstrap_statistics_pairwise,
                                Rcpp::Named("probs") = 1 - level));
  }
  if (progress) {
    REprintf("\nComputing confidence intervals...\n");
  }
  if (interval) {
    std::vector<double> lower(m);
    std::vector<double> upper(m);
    for (int i = 0; i < m; ++i) {
      Rcpp::checkUserInterrupt();
      Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(1, x.cols());
      lhs(pairs[i][0]) = 1;
      lhs(pairs[i][1]) = -1;
      std::array<double, 2> ci =
        pair_confidence_interval_gbd(theta_hat, x, c, lhs,
                                     threshold, estimate[i], cutoff);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }
    result["lower"] = lower;
    result["upper"] = upper;
  }
  result["p.adj"] = adj_pvalues;
  result["k"] = k;
  result["level"] = level;
  result["method"] = method;
  result["B"] = bootstrap_statistics_pairwise.size();
  result.attr("class") = Rcpp::CharacterVector({"pairwise", "elgbd"});
  return result;
}
