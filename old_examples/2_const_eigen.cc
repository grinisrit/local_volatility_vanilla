#include <iostream>
#include <vector>

#define EIGEN_NO_AUTOMATIC_RESIZING 1
#define EIGEN_DONT_ALIGN 1
#define EIGEN_NO_DEBUG 1
#define EIGEN_UNROLLING_LIMIT 0
#define EIGEN_DONT_VECTORIZE 1


#include <Eigen/Dense>


constexpr size_t IN = 4, OUT = 4;
using Matrix = Eigen::Matrix<double, IN, OUT>; //Eigen::MatrixXd;
using DynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

double prod_norm(const Matrix& X, const Matrix& Y) {
  return (X*Y).sum();
}

void __enzyme_autodiff(void*, int, const Matrix&, const Matrix&, int, const Matrix&);
int enzyme_const, enzyme_dup;


auto main() -> int {
  double h = 0.001;
  int row = 3, col = 2;

  DynamicMatrix A,B; A.resize(4,4); B.resize(4,2); 
  B << 1,2,3,4,5,6,7,8;

  std::cout << "A:\n" << A << std::endl;
  std::cout << "B:\n" << B << std::endl;

  std::vector<double> inputs(A.rows() * A.cols());
  for (int i=0; i < inputs.size(); i++) inputs[i] = i*i;

  double w = 2.;
  double m = 1.;
  const Matrix W = Matrix::Constant(IN, OUT, w);
  const Matrix M = Matrix::Constant(IN, OUT, m);
  const Matrix Wp = Matrix::Constant(IN, OUT, 0.0);
  
  __enzyme_autodiff((void*)prod_norm, 
    enzyme_dup, W, Wp, 
    enzyme_const, M);


  std::cout << Wp << std::endl;

  return 0;
}