/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Iinclude -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp
*/

#include <random>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "test_simplex.h"
#include <ext/uniform_inside_tetrahedron_distribution.h>

int
main()
{
  using array_t = __gnu_test::uniform_inside_tetrahedron_distribution<>::result_type;
  using real_t = array_t::value_type;

  std::cout.precision(std::numeric_limits<real_t>::digits10);
  auto w = 6 + std::cout.precision();

  __gnu_test::uniform_inside_tetrahedron_distribution<> uit;

  std::cout << "uit = " << uit << '\n';

  tetra<real_t> tetra0;
  tetra<real_t> tetra1({6, 6, 6}, {-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

  std::random_device rd;
  std::mt19937 gen(rd());

  for (int i = 0; i < 100000; ++i)
    {
      auto bary = uit(gen);
      std::cout << ' ' << std::setw(w) << bary[0]
		<< ' ' << std::setw(w) << bary[1]
		<< ' ' << std::setw(w) << bary[2]
		<< ' ' << std::setw(w) << bary[3] << '\n';
    }

  std::ofstream iv("tetras.iv");
  iv << "#Inventor V2.1 ascii\n";
  iv << "Separator {\n";

  iv << "  BaseColor {\n";
  iv << "    rgb 1.0 0.0 0.0\n";
  iv << "  }\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";
  for (int i = 0; i < 100000; ++i)
    {
      auto bary = uit(gen);
      iv << tetra0(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";

  iv << "  BaseColor {\n";
  iv << "    rgb 0.0 0.0 1.0\n";
  iv << "  }\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";
  for (int i = 0; i < 100000; ++i)
    {
      auto bary = uit(gen);
      iv << tetra1(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";

  iv << "}\n";
}

