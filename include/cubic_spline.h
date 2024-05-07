#pragma once

#include <vector>


inline size_t search_sorted(const std::vector<double>& a, const double b) {
  const size_t n = a.size();
  size_t idx = 0;
  size_t pa = 0;
  size_t pb = 0;

  while (pb < 1){
    if (pa < n-1 && a[pa] < b){
      pa += 1;
    } else {
      idx = pa;
      pb += 1;
    }
  }
  return idx;
} 


class CubicSpline {
  public:
  CubicSpline(const std::vector<double>& x, const std::vector<double>& y) {
    const size_t n = x.size() - 1;

    std::vector<double> h;
    h.reserve(n);
    for(size_t i = 0; i < n; i++) {
      h.push_back(x[i+1] - x[i]);
    }

    std::vector<double> alpha;
    alpha.reserve(n-1);
    for(size_t i = 0; i < n-1; i++) {
      alpha.push_back(
        3 * ((y[i+2] - y[i+1]) / h[i+1] - (y[i+1] - y[i]) / h[i])
      );
    }

    std::vector<double> c = std::vector<double>(n+1, 0.);
    std::vector<double> ell = std::vector<double>(n+1, 1.);
    std::vector<double> mu = std::vector<double>(n, 0.);
    std::vector<double> z = std::vector<double>(n+1, 0.);

    for(size_t i = 1; i < n; i++) {
      ell[i] = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1];
      mu[i] = h[i] / ell[i];
      z[i] = (alpha[i-1] - h[i-1] * z[i-1]) / ell[i];
    }

    for(size_t i = 0; i < n; i++) {
      size_t j = n - 1 - i;
      c[j] = z[j] - mu[j] * c[j+1];
    }

    x0_.reserve(n);
    for(size_t i = 0; i < n; i++) {
      x0_.push_back(x[i+1]);
    }

    a_.reserve(n);
    for(size_t i = 0; i < n; i++) {
      a_.push_back(y[i+1]);
    }

    b_.reserve(n);
    for(size_t i = 0; i < n; i++) {
      b_.push_back(
        (y[i+1] - y[i]) / h[i] + (c[i] + 2 * c[i+1]) * h[i] / 3
      );
    }

    c_.reserve(n);
    for(size_t i = 0; i < n; i++) {
      c_.push_back(c[i+1]);
    }

    d_.reserve(n);
    for(size_t i = 0; i < n; i++) {
      d_.push_back(
        (c[i+1] - c[i]) / (3*h[i])
      );
    }    
  }

  double apply(double x) const {
    size_t ix = search_sorted(x0_, x);
    double dx = x - x0_[ix];
    return a_[ix] + (b_[ix] + (c_[ix] + d_[ix] * dx) * dx) * dx;
  }

  private:
  std::vector<double> x0_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
  std::vector<double> d_;

};