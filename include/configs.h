#pragma once

#include <vector>

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