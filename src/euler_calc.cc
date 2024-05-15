#include "euler_calc.h"

#include "spline.h"

#include "deviates.h"


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
