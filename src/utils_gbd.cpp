#include "utils_gbd.h"

Eigen::MatrixXd g_gbd(const Eigen::Ref<const Eigen::VectorXd>& theta,
                      const Eigen::Ref<const Eigen::MatrixXd>& x,
                      const Eigen::Ref<const Eigen::MatrixXd>& c) {
  return x - (c.array().rowwise() * theta.array().transpose()).matrix();
}

Eigen::MatrixXd cov_gbd(const Eigen::Ref<const Eigen::MatrixXd>& x,
                        const Eigen::Ref<const Eigen::MatrixXd>& c) {
  // Estimating function
  Eigen::MatrixXd g =
    g_gbd(x.array().colwise().sum() / c.array().colwise().sum(), x, c);
  // Covariance estimate
  return (g.transpose() * g) / x.rows();
}

Eigen::VectorXd lambda2theta_gbd(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const double gamma) {
  Eigen::VectorXd ngradient =
    (PSEUDO_LOG::dp(Eigen::VectorXd::Ones(g.rows()) + g * lambda).matrix().asDiagonal() * c)
  .array().colwise().sum().transpose() * lambda.array();
  return theta + gamma * ngradient;
}

void lambda2theta_void(
    const Eigen::Ref<const Eigen::VectorXd>& lambda,
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& g,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const double gamma) {
  theta +=
    gamma * ((PSEUDO_LOG::dp(
    Eigen::VectorXd::Ones(
      g.rows()) + g * lambda).matrix().asDiagonal() * c)
               .array().colwise().sum().transpose() *
      lambda.array()).matrix();
}

Eigen::VectorXd approx_lambda_gbd(
    const Eigen::Ref<const Eigen::MatrixXd>& g0,
    const Eigen::Ref<const Eigen::MatrixXd>& c,
    const Eigen::Ref<const Eigen::VectorXd>& theta0,
    const Eigen::Ref<const Eigen::VectorXd>& theta1,
    const Eigen::Ref<const Eigen::VectorXd>& lambda0) {
  Eigen::ArrayXd&& arg = Eigen::VectorXd::Ones(g0.rows()) + g0 * lambda0;
  Eigen::ArrayXd&& denominator = Eigen::pow(arg, 2);

  // LHS
  Eigen::MatrixXd&& LHS =
    g0.transpose() * (g0.array().colwise() / denominator).matrix();

  // RHS
  const Eigen::MatrixXd I_RHS =
    ((c.array().colwise() / arg).colwise().sum()).matrix().asDiagonal();
  const Eigen::MatrixXd J_RHS =
    (g0.array().colwise() / denominator).matrix().transpose() *
    (c.array().rowwise() * lambda0.array().transpose()).matrix();
  Eigen::MatrixXd&& RHS = -I_RHS + J_RHS;

  // Jacobian matrix
  Eigen::MatrixXd&& jacobian = LHS.ldlt().solve(RHS);

  // Linear approximation for lambda1
  return lambda0 + jacobian * (theta1 - theta0);
}

Eigen::MatrixXd rmvn(const Eigen::MatrixXd& x, const int n) {
  // Generate standard multivariate normal random vectors (n by p matrix)
  Eigen::MatrixXd I(n, x.cols());
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < n; ++i) {
      I(i, j) = R::rnorm(0, 1.0);
    }
  }
  // Get the square root matrix of the covariance matrix
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(x);
  // Return the target normal random vectors(n by p matrix)
  return I * es.operatorSqrt();
}

minEL test_gbd_EL(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                  const Eigen::Ref<const Eigen::MatrixXd>& x,
                  const Eigen::Ref<const Eigen::MatrixXd>& c,
                  const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                  const Eigen::Ref<const Eigen::VectorXd>& rhs,
                  const int maxit,
                  const double abstol,
                  const double threshold) {
  /// Initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(theta0, lhs, rhs);
  // Estimating function
  Eigen::MatrixXd g = g_gbd(theta, x, c);
  // Evaluation
  Eigen::VectorXd lambda = EL(g, maxit, abstol, threshold).lambda;
  // For current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// Minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // Proposed value for theta
  while (!convergence && iterations != maxit) {
    // Update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // Update g
    Eigen::MatrixXd g_tmp = g_gbd(theta_tmp, x, c);
    // Update lambda
    EL eval(g_tmp, maxit, abstol, threshold);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      return {theta, lambda, f1, iterations, convergence};
    }

    // Update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    // Step halving to ensure that the updated function value be
    // trictly less than the current function value.
    while (f0 < f1) {
      // Reduce step size
      gamma /= 2;
      // Propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // Propose new lambda
      g_tmp = g_gbd(theta_tmp, x, c);
      lambda_tmp = EL(g_tmp, maxit, abstol, threshold).lambda;
      if (gamma < abstol) {
        return {theta, lambda, f0, iterations, convergence};
      }
      // Propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // Update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // Convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return {theta, lambda, f1, iterations, convergence};
}

double test_nlogLR(const Eigen::Ref<const Eigen::VectorXd>& theta0,
                   const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol,
                   const double threshold) {
  /// Initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(theta0, lhs, rhs);
  // Estimating function
  Eigen::MatrixXd g = g_gbd(theta, x, c);
  // Evaluation
  Eigen::VectorXd lambda = EL(g, maxit, abstol, threshold).lambda;
  // For current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// Minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // Proposed value for theta
  while (!convergence && iterations != maxit) {
    // Update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // Update g
    Eigen::MatrixXd g_tmp = g_gbd(theta_tmp, x, c);
    // Update lambda
    EL eval(g_tmp, maxit, abstol, threshold);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      return f1;
    }

    // Update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);

    // Step halving to ensure that the updated function value be
    // strictly less than the current function value.
    while (f0 < f1) {
      // Reduce step size
      gamma /= 2;
      // Propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // Propose new lambda
      g_tmp = g_gbd(theta_tmp, x, c);
      lambda_tmp = EL(g_tmp, maxit, abstol, threshold).lambda;
      if (gamma < abstol) {
        return f0;
      }
      // Propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // Update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // Convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return f1;
}

double test_nlogLR(const Eigen::Ref<const Eigen::MatrixXd>& x,
                   const Eigen::Ref<const Eigen::MatrixXd>& c,
                   const Eigen::Ref<const Eigen::MatrixXd>& lhs,
                   const Eigen::Ref<const Eigen::VectorXd>& rhs,
                   const int maxit,
                   const double abstol,
                   const double threshold) {
  /// Initialization ///
  // Constraint imposed on the initial value by projection.
  // The initial value is given as treatment means.
  Eigen::VectorXd theta =
    linear_projection(x.array().colwise().sum() / c.array().colwise().sum(),
                      lhs, rhs);
  // Estimating function
  Eigen::MatrixXd g = g_gbd(theta, x, c);
  // Evaluation
  Eigen::VectorXd lambda = EL(g, maxit, abstol, threshold).lambda;
  // For current function value(-logLR)
  double f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g.rows()) + g * lambda);

  /// Minimization(projected gradient descent) ///
  double gamma = 1.0 / (c.colwise().sum().mean());    // step size
  bool convergence = false;
  int iterations = 0;
  // Proposed value for theta
  while (!convergence && iterations != maxit) {
    // Update parameter by GD with lambda fixed -> projection
    Eigen::VectorXd theta_tmp = theta;
    lambda2theta_void(lambda, theta_tmp, g, c, gamma);
    linear_projection_void(theta_tmp, lhs, rhs);
    // Update g
    Eigen::MatrixXd g_tmp = g_gbd(theta_tmp, x, c);
    // Update lambda
    EL eval(g_tmp, maxit, abstol, threshold);
    Eigen::VectorXd lambda_tmp = eval.lambda;
    if (!eval.convergence && iterations > 9) {
      return f1;
    }
    // Update function value
    double f0 = f1;
    f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    // Step halving to ensure that the updated function value be
    // strictly less than the current function value.
    while (f0 < f1) {
      // Reduce step size
      gamma /= 2;
      // Propose new theta
      theta_tmp = theta;
      lambda2theta_void(lambda, theta_tmp, g, c, gamma);
      linear_projection_void(theta_tmp, lhs, rhs);
      // Propose new lambda
      g_tmp = g_gbd(theta_tmp, x, c);
      lambda_tmp = EL(g_tmp, maxit, abstol, threshold).lambda;
      if (gamma < abstol) {
        return f0;
      }
      // Propose new function value
      f1 = PSEUDO_LOG::sum(Eigen::VectorXd::Ones(g_tmp.rows()) + g_tmp * lambda_tmp);
    }

    // Update parameters
    theta = std::move(theta_tmp);
    lambda = std::move(lambda_tmp);
    g = std::move(g_tmp);
    ++iterations;

    // Convergence check
    if (f0 - f1 < abstol) {
      convergence = true;
    }
  }

  return f1;
}
