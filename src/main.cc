#include <iostream>
#include <vector>

#include "rapidcsv.h"
#include "deviates.h"

#include "configs.h"
#include "vol_surface.h"

#include "spline.h"

inline const double DT = 2.7E-3;

extern int enzyme_const, enzyme_dup;
      
template <typename Retval, typename... Args>
Retval __enzyme_autodiff(Retval (*)(Args...), auto...);

double calc_pv(const VolMatrix& sigmas, const MarketDataConfig& market_config, const VanillaCallTradeConfig& trade_config) { 

  tk::mat_double strikes;
  strikes.reserve(1);
  strikes.push_back(market_config.strikes);

  const tk::spline strike_splines = tk::spline{
    strikes, sigmas
    };

  const VolSurface vol_surface = VolSurface(market_config, sigmas);
  
  Normaldev normal_dev{0., 1., 10};

  double pv = 0.;
  for (size_t n = 0; n < trade_config.N; n++) {
    double S = vol_surface.flat_forward();
    double T = 0.;
    for(size_t ttm = 0; ttm < trade_config.TTM; ttm++) {
      T = T + DT;
      S = S + vol_surface.get_local_vol(S,T, strike_splines) * S * sqrt(DT) * normal_dev.dev();
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
  /**/
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