#define EIGEN_NO_AUTOMATIC_RESIZING 1
#define EIGEN_DONT_ALIGN 1
#define EIGEN_NO_DEBUG 1
#define EIGEN_UNROLLING_LIMIT 0
#define EIGEN_DONT_VECTORIZE 1


#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;

extern "C" {
    extern double __enzyme_autodiff(void*, const MatrixXd* __restrict W, const MatrixXd* __restrict Wp, const MatrixXd* __restrict M, const MatrixXd* __restrict Mp);
}

__attribute__((noinline))
static double matvec(const MatrixXd* __restrict W, const MatrixXd* __restrict M) {
  MatrixXd diff = *W-*M;
  return (diff*diff).sum();
}

auto main() -> int {
  size_t IN = 40, OUT = 30, NUM = 50;
  MatrixXd W = Eigen::MatrixXd::Constant(IN, OUT, 1.0);
  MatrixXd M = Eigen::MatrixXd::Constant(IN, OUT, 2.0);
  
  MatrixXd Wp = Eigen::MatrixXd::Constant(IN, OUT, 0.0);
  MatrixXd Mp = Eigen::MatrixXd::Constant(IN, OUT, 0.0);

  __enzyme_autodiff((void*)matvec, &W, &Wp, &M, &Mp);

  return 0;
}