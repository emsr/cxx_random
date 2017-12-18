#include <array>
#include <cmath>
#include <iostream>
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

  for (int i = 0; i < 1.0e8; ++i)
  {
    const int number = distribution(engine);
    hist << number;
  }
  int k = 1;
  for (const auto& bin : hist)
    std::cout << hist.lower_bound(k++) << ' ' << bin << ' ' << bin / 1.0e8 << '\n';
}
