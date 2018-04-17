/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp
*/

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "test_simplex.h"
#include "uniform_inside_triangle_distribution.h"


int
main()
{
  using array_t = __gnu_test::uniform_inside_triangle_distribution<>::result_type;
  using real_t = array_t::value_type;

  std::cout.precision(std::numeric_limits<real_t>::digits10);
  auto w = 6 + std::cout.precision();

  __gnu_test::uniform_inside_triangle_distribution<> uit;

  std::cout << "uit = " << uit << '\n';

  tri<real_t> tri0;
  tri<real_t> tri1({-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

  std::random_device rd;
  std::mt19937 gen(rd());

  for (int i = 0; i < 1000; ++i)
    {
      auto bary = uit(gen);
      std::cout << ' ' << std::setw(w) << bary[0]
		<< ' ' << std::setw(w) << bary[1]
		<< ' ' << std::setw(w) << bary[2] << '\n';
    }

  std::ofstream iv("tris.iv");
  iv << "#Inventor V2.1 ascii\n";
  iv << "Separator {\n";

  iv << "  BaseColor {\n";
  iv << "    rgb 1.0 0.0 0.0\n";
  iv << "  }\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";
  for (int i = 0; i < 1000; ++i)
    {
      auto bary = uit(gen);
      iv << tri0(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";

  iv << "  BaseColor {\n";
  iv << "    rgb 0.0 0.0 1.0\n";
  iv << "  }\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";
  for (int i = 0; i < 1000; ++i)
    {
      auto bary = uit(gen);
      iv << tri1(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";

  iv << "}\n";
}

