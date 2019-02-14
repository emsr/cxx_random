/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Iinclude -Wall -Wextra -o test_dirichlet_distribution test_dirichlet_distribution.cpp
*/

#include <ext/random>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <ext/dirichlet_distribution.h>
#include <ext/simplex.h>
#include "test_simplex.h"

template<typename RealType = double>
  void
  test_dirichlet_triangle()
  {
    std::cout.precision(std::numeric_limits<RealType>::digits10);
    auto w = 6 + std::cout.precision();

    __gnu_cxx::dirichlet_distribution<3> dd{1.0, 1.0, 3.0};

    std::cout << "dd = " << dd << '\n';

    tri<RealType> tri0;
    tri<RealType> tri1({-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < 1000; ++i)
      {
	auto bary = dd(gen);
	std::cout << ' ' << std::setw(w) << bary[0]
		  << ' ' << std::setw(w) << bary[1]
		  << ' ' << std::setw(w) << bary[2] << '\n';
      }

    std::ofstream iv("dirichlet_tris.iv");
    iv << "#Inventor V2.1 ascii\n";
    iv << "Separator {\n";

    iv << "  BaseColor {\n";
    iv << "    rgb 1.0 0.0 0.0\n";
    iv << "  }\n";
    iv << "  Coordinate3 {\n";
    iv << "    point [\n";
    for (int i = 0; i < 1000; ++i)
      {
	auto bary = dd(gen);
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
	auto bary = dd(gen);
	iv << tri1(bary) << ",\n";
      }
    iv << "    ]\n";
    iv << "  }\n";
    iv << "  PointSet {}\n";

    iv << "}\n";
  }

template<typename RealType = double>
  void
  test_dirichlet_tetrahedron()
  {
    std::cout.precision(std::numeric_limits<RealType>::digits10);
    auto w = 6 + std::cout.precision();

    __gnu_cxx::dirichlet_distribution<4> dd{5.0, 1.0, 1.0, 1.0};

    std::cout << "dd = " << dd << '\n';

    tetra<RealType> tetra0;
    tetra<RealType> tetra1({6, 6, 6}, {-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < 100000; ++i)
      {
	auto bary = dd(gen);
	std::cout << ' ' << std::setw(w) << bary[0]
		  << ' ' << std::setw(w) << bary[1]
		  << ' ' << std::setw(w) << bary[2] << '\n';
      }

    std::ofstream iv("dirichlet_tetras.iv");
    iv << "#Inventor V2.1 ascii\n";
    iv << "Separator {\n";

    iv << "  BaseColor {\n";
    iv << "    rgb 1.0 0.0 0.0\n";
    iv << "  }\n";
    iv << "  Coordinate3 {\n";
    iv << "    point [\n";
    for (int i = 0; i < 100000; ++i)
      {
	auto bary = dd(gen);
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
	auto bary = dd(gen);
	iv << tetra1(bary) << ",\n";
      }
    iv << "    ]\n";
    iv << "  }\n";
    iv << "  PointSet {}\n";

    iv << "}\n";
  }

int
main()
{
  test_dirichlet_triangle();
  test_dirichlet_tetrahedron();
}
