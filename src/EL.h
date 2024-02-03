#ifndef EL_H_
#define EL_H_

#include "eigen_config.h"
#include <RcppEigen.h>
#include <cmath>
#include "utils.h"

struct minEL {
  Eigen::VectorXd par;
  Eigen::VectorXd lambda;
  double nlogLR;
  int iterations;
  bool convergence;
};

class EL
{
private:
  const int n;
public:
  Eigen::VectorXd lambda;
  double nlogLR = 0;
  int iterations = 1;
  bool convergence = false;

  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const int maxit,
     const double abstol,
     const double threshold);

  EL(const Eigen::Ref<const Eigen::MatrixXd>& g,
     const Eigen::Ref<const Eigen::ArrayXd>& w,
     const int maxit,
     const double abstol,
     const double threshold);
};

class PSEUDO_LOG
{
public:
  Eigen::ArrayXd dplog;
  Eigen::ArrayXd sqrt_neg_d2plog;
  double plog_sum = 0;
  PSEUDO_LOG(Eigen::VectorXd&& x);
  static Eigen::ArrayXd plog(Eigen::ArrayXd&& x);
  static double sum(Eigen::VectorXd&& x);
  static Eigen::ArrayXd dp(Eigen::VectorXd&& x);
};
#endif
