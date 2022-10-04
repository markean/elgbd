#include "utils.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int get_max_threads() {
  #ifdef _OPENMP
  return omp_get_max_threads();
  #else
  return 1;
  #endif
}

Eigen::VectorXd linear_projection(
    const Eigen::Ref<const Eigen::VectorXd>& theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs) {
  return theta -
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * (lhs * theta - rhs);
}
void linear_projection_void(
    Eigen::Ref<Eigen::VectorXd> theta,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::VectorXd>& rhs) {
  theta -=
    lhs.transpose() * (lhs * lhs.transpose()).inverse() * (lhs * theta - rhs);
}

Eigen::MatrixXd bootstrap_sample(
    const Eigen::Ref<const Eigen::MatrixXd>& x,
    const Eigen::Ref<const Eigen::ArrayXi>& index) {
  Eigen::MatrixXd out(x.rows(), x.cols());
  for (int i = 0; i < x.rows(); ++i) {
    out.row(i) = x.row(index(i));
  }
  return out;
}
