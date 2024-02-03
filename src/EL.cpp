#include "EL.h"

// Constructor for EL class (evaluation)
EL::EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
       const int maxit,
       const double abstol,
       const double threshold)
  : n{static_cast<int>(g.rows())}
{
  // Maximization
  lambda = (g.transpose() * g).ldlt().solve(g.colwise().sum());
  while (!convergence && iterations != maxit) {
    // Plog class
    PSEUDO_LOG log_tmp(Eigen::VectorXd::Ones(n) + g * lambda);
    // J matrix
    const Eigen::MatrixXd J = g.array().colwise() * log_tmp.sqrt_neg_d2plog;
    // Propose new lambda by NR method with least square
    Eigen::VectorXd step =
      (J.transpose() * J).ldlt().solve(
          J.transpose() * (log_tmp.dplog / log_tmp.sqrt_neg_d2plog).matrix());
    // Update function value
    nlogLR = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(n) + g * (lambda + step));
    // Step halving to ensure increase in function value
    double gamma = 1.0;
    while (nlogLR < log_tmp.plog_sum) {
      gamma /= 2;
      if (gamma < abstol) {
        break;
      }
      nlogLR =
        PSEUDO_LOG::sum(
          Eigen::VectorXd::Ones(n) + g * (lambda + gamma * step));
    }
    /* If the step halving is not successful (possibly due to the convex
     * hull constraint), terminate the maximization with the current values
     * without further updates.
     */
    if (gamma < abstol) {
      nlogLR = log_tmp.plog_sum;
      break;
    }
    // Otherwise, update lambda and check for convergence
    lambda += gamma * step;
    if (nlogLR - log_tmp.plog_sum < abstol) {
      convergence = true;
    } else {
      ++iterations;
    }
  }
}

// Constructor for PSEUDO_LOG class
PSEUDO_LOG::PSEUDO_LOG(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;

  dplog.resize(x.size());
  sqrt_neg_d2plog.resize(x.size());

  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      dplog[i] = a2 + 2.0 * a3 * x[i];
      sqrt_neg_d2plog[i] = a2 / 2.0;
      plog_sum += a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      dplog[i] = 1.0 / x[i];
      sqrt_neg_d2plog[i] = 1.0 / x[i];
      plog_sum += log(x[i]);
    }
  }
}

// Pseudo log function
Eigen::ArrayXd PSEUDO_LOG::plog(Eigen::ArrayXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;
  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      x[i] = a1 + a2 * x[i] + a3 * x[i] * x[i];
    } else {
      x[i] = log(x[i]);
    }
  }
  return x;
}

double PSEUDO_LOG::sum(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = -log(n) - 1.5;
  static const double a2 = 2.0 * n;
  static const double a3 = -0.5 * n * n;
  double out = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    out += x[i] < a0 ? a1 + a2 * x[i] + a3 * x[i] * x[i] : log(x[i]);
  }
  return out;
}

Eigen::ArrayXd PSEUDO_LOG::dp(Eigen::VectorXd&& x) {
  static const double n = static_cast<double>(x.size());
  static const double a0 = 1.0 / n;
  static const double a1 = 2.0 * n;
  static const double a2 = -1.0 * n * n;
  for (unsigned int i = 0; i < x.size(); ++i) {
    if (x[i] < a0) {
      x[i] = a1 + a2 * x[i];
    } else {
      x[i] = 1.0 / x[i];
    }
  }
  return x;
}
