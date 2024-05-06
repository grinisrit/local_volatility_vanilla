#include <iostream>
#include <vector>
#include <omp.h>

double sum_vec(const std::vector<double>& vec) {
  double res = 0.;
  const size_t n = vec.size();
  for (size_t i = 0; i<n; i++) {
    res += vec.at(i);
  }
  return res;
}

double square_norm(const std::vector<double>& vec) {
  const size_t n = vec.size();
  std::vector<double> res =  std::vector(n,0.);

#pragma omp parallel for
  for (size_t i = 0; i<n; i++) {
    res[i] = vec.at(i) * vec.at(i);
  }
  return sum_vec(res);
}

void __enzyme_autodiff(void*, const std::vector<double>&, const std::vector<double>&);


auto main(int argc, char *argv[]) -> int {
  int n = 1000;
  double x = 20.;
  if (argc > 1) {
    n = std::atoi(argv[1]);
    if (argc > 2) {
      x = std::atof(argv[2]);
    }
  }

  const std::vector<double> vec = std::vector(n,x);
  const std::vector<double> grad_vec = std::vector(n, 0.);
  std::cout << "Num of threads: " << omp_get_max_threads() << std::endl;
  std::cout << "Square norm: " << square_norm(vec) << std::endl;

  __enzyme_autodiff((void*)square_norm,  vec, grad_vec);

  std::cout << " Gradient sum = " << sum_vec(grad_vec) << std::endl;

  return 0;
}