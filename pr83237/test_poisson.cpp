#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include "../histogram.h"

int
main()
{
  // The problem turned out to be independent on the engine
  std::mt19937_64 engine;

  // Set fixed seed for easy reproducibility
  // The problem turned out to be independent on seed
  engine.seed(1);
  std::poisson_distribution<int> distribution(157.17);
  histogram<int> hist(140, 90, 230);

  const auto max = 1.0e8;
  for (int i = 0; i < max; ++i)
  {
    const int number = distribution(engine);
    hist << number;
  }
  int k = 0;
  for (const auto& bin : hist)
    std::cout << ' ' << std::setw(4) << hist.lower_bound(++k)
	      << ' ' << std::setw(8) << bin
	      << ' ' << std::setw(12) << bin / max
	      << '\n';
}
