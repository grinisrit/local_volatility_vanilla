#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "rapidcsv.h"
#include "deviates.h"


inline const double EPS = 1E-3;
inline const double TOL = 1E-9;
inline const double DT = 2.7E-3;

using VolMatrix = std::vector<std::vector<double>>;

extern int enzyme_const, enzyme_dup;
      
template <typename Retval, typename... Args>
Retval __enzyme_autodiff(Retval (*)(Args...), auto...);

inline double normal_cdf(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}


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


class CubicSpline1D {
  public:
  CubicSpline1D(const std::vector<double>& x, const std::vector<double>& y) {
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


struct VanillaCallTradeConfig {
  const double K;
  const size_t N;
  const size_t TTM;
};

struct MarketDataConfig
{
  const std::vector<double> ttms;
  const std::vector<double> strikes;
  const double fwd; 
};


class VolSurface {

  double clip_ttm(double ttm) const {
    ttm = std::max(ttm, ttms_[0]); // assuming first quotes are EOD current date 
    ttm = std::min(ttm, ttms_.back());
    return ttm;
  }

  double clip_F(double F) const {
    F = std::max(F, min_K_); 
    F = std::min(F, max_K_);
    return F;
  }

  std::tuple<size_t, size_t> get_ttm_interval(double ttm) const {
    ttm = std::min(ttm, ttms_.back());
    size_t t1 = 0; 
    size_t t2 = 1;

    for(size_t i = 1; i < n_ttms_; i++) {
      const auto& ttm_ = ttms_[i];
      if (ttm_ > ttm) {
        t1 = i-1;
        t2 = i;
        break;
      }
    }
    return std::make_tuple(t1, t2);
  }

  double compute_sigma(double ttm, size_t t1, size_t t2, double F) const {
    const double sigma1 = std::max(strike_splines_[t1].apply(F), TOL); // const fwd yield
    const double sigma2 = std::max(strike_splines_[t2].apply(F), TOL);
    const double ttm1 = ttms_[t1]; 
    const double ttm2 = ttms_[t2];
    const double dttm = ttm2 - ttm1;


    const double tot_var = (ttm2 * sigma2 * sigma2) * (ttm2 - ttm) / dttm + 
      (ttm1 * sigma1 * sigma1) * (ttm - ttm1) / dttm;

    return sqrt(tot_var / ttm);
  }

  double bsm_pv(double ttm, double K, double sigma) const {
    const double d1 = (log(fwd_ / K) + (sigma*sigma / 2) * ttm) / (sigma * sqrt(ttm));
    const double d2 = d1 - sigma * sqrt(ttm);
    return fwd_ * normal_cdf(d1) - K * normal_cdf(d2);
  }

  double compute_C(double ttm, double F) const {
    const auto& [t1, t2] = get_ttm_interval(ttm);
    const double sigma = compute_sigma(ttm, t1, t2, F);
    return bsm_pv(ttm, F, sigma);
  }

  double finite_diff_dupire(double ttm, double F) const {
    const double C = compute_C(ttm, F);
    const double C_Kl = compute_C(ttm, clip_F(F - EPS));
    const double C_Ku = compute_C(ttm, clip_F(F + EPS));
    const double d2C = C_Ku - 2 * C + C_Kl;
    const double d2C_dK2 = std::max(d2C / (EPS*EPS), TOL); // butterly arbitrage

    const double C_tl = compute_C(clip_ttm(ttm - EPS), F);
    const double C_tu = compute_C(clip_ttm(ttm + EPS), F);
    const double dC = C_tu - C_tl;
    const double dC_dT = std::max(dC / EPS, TOL); // calendar arbitrage 

    const double loc_vol2 = (2*dC_dT) / (F*F*d2C_dK2);
    
    return sqrt(loc_vol2);
  }

  public:
  VolSurface(
    const MarketDataConfig& market_config, 
    const VolMatrix& sigmas) : fwd_(market_config.fwd), n_ttms_(market_config.ttms.size()) {
  
  ttms_.reserve(n_ttms_);
  for(const auto& ttm: market_config.ttms) {
    ttms_.push_back(ttm);
  }

  strike_splines_.reserve(n_ttms_);
  for(size_t i = 0; i < n_ttms_; i++) {
    strike_splines_.emplace_back(market_config.strikes, sigmas[i]);
  }

  max_K_ = market_config.strikes.back();
  min_K_ = market_config.strikes.front();
  }

  double get_local_vol(double K, double ttm) const {
    ttm = clip_ttm(ttm);
    K = clip_F(K);
    return finite_diff_dupire(ttm, K);
  }

  double flat_forward() const {
    return fwd_;
  }

  private:
  double fwd_;
  double max_K_;
  double min_K_;
  size_t n_ttms_;
  std::vector<CubicSpline1D> strike_splines_;
  std::vector<double> ttms_;

};

double calc_pv(const VolMatrix& sigmas, const MarketDataConfig& market_config, const VanillaCallTradeConfig& trade_config) {  

  const VolSurface vol_surface = VolSurface(market_config, sigmas);
  
  Normaldev normal_dev{0., 1., 10};

  double pv = 0.;
  for (size_t n = 0; n < trade_config.N; n++) {
    double S = vol_surface.flat_forward();
    double T = 0.;
    for(size_t ttm = 0; ttm < trade_config.TTM; ttm++) {
      T = T + DT;
      S = S + vol_surface.get_local_vol(S,T) * S * sqrt(DT) * normal_dev.dev();
    }
    double payoff = std::max(S - trade_config.K, 0.);
    pv = pv + payoff;
  }
  return pv/trade_config.N;              
}


auto main(int argc, char *argv[]) -> int {
  size_t N_PATHS = 1000;
  size_t N_DAYS = 300;
  if (argc > 1) {
    N_PATHS = (size_t) std::atoi(argv[1]);
    if (argc > 2) {
      N_DAYS = std::min((size_t) std::atoi(argv[2]), N_DAYS);
    }
  }

  rapidcsv::Document fwd_csv("../fwd.csv", rapidcsv::LabelParams(-1, -1));
  std::vector<double> ttms = fwd_csv.GetColumn<double>(0);
  const size_t n_ttms = ttms.size();

  double fwd = fwd_csv.GetColumn<double>(1).at(0);

  rapidcsv::Document impl_vol_csv("../impl_vol.csv", rapidcsv::LabelParams(-1, -1));

  std::vector<double> strikes = impl_vol_csv.GetColumn<double>(0);
  const size_t n_strikes = strikes.size();
  
  VolMatrix sigmas;
  sigmas.reserve(n_ttms);
  for (size_t col = 1; col < n_ttms + 1; col++) {
    std::vector<double> smile = impl_vol_csv.GetColumn<double>(col);
    sigmas.push_back(smile);
  }

  const MarketDataConfig market_config = MarketDataConfig{ttms, strikes, fwd};
  
  double K = 1.1 * fwd;
  const VanillaCallTradeConfig trade_config = VanillaCallTradeConfig{K, N_PATHS, N_DAYS};

  double pv = calc_pv(sigmas, market_config, trade_config);
  std::cout << "PV: " << pv << std::endl;

  VolMatrix vegas;
  vegas.reserve(n_ttms);
  for (size_t i = 0; i < n_ttms; i++) {
    vegas.push_back(std::vector<double>(n_strikes, 0.));
  }

 
  __enzyme_autodiff(calc_pv, 
      enzyme_dup, &sigmas, &vegas,
      enzyme_const, &market_config, 
      enzyme_const, &trade_config);

  
  std::cout << "Vegas:\n" << 
    vegas[0][0] << " " <<  vegas[0][1] <<  "...\n" <<
    vegas[1][0] << " " <<  vegas[1][1] <<  "..." << std::endl;

  std::ofstream vegas_csv("vegas.csv");
  for(const auto& vega_t: vegas) {
    for(const auto& vega: vega_t) {
      vegas_csv << vega << ",";
    }
    vegas_csv << "\n";
  }
  vegas_csv.close();

}