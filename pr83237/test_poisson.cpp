#include <array>
#include <cmath>
#include <iostream>
#include <random>

int
main()
{
  // The problem turned out to be independent on the engine
  std::mt19937_64 engine;

  // Set fixed seed for easy reproducibility
  // The problem turned out to be independent on seed
  engine.seed(1);
  std::poisson_distribution<int> distribution(157.17);

  for (int i = 0; i < 1E8; i++) {
    const int number = distribution(engine);
    std::cout << number << '\n';
  }
}