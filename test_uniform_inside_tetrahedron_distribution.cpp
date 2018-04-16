/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp
*/

#include <random>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "simplex.h"

template<class RealTp = double>
  constexpr std::array<std::array<RealTp, 3>, 4>
  tetra0{{{1, 1, 1}, {0, 0, 1}, {1, 0, 0}, {0, 1, 0}}};

template<class RealTp = double>
  constexpr std::array<std::array<RealTp, 3>, 4>
  tetra1{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

template<class RealTp = double>
  constexpr std::array<std::array<RealTp, 3>, 4>
  tetra2{{{1, 0, 1}, {0, 0, 1}, {1, 0, 0}, {1, 1, 1}}};

template<class RealTp = double>
  constexpr std::array<std::array<RealTp, 3>, 4>
  tetra3{{{0, 1, 1}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}};

template<class RealTp = double>
  constexpr std::array<std::array<RealTp, 3>, 4>
  tetra4{{{1, 1, 0}, {0, 1, 0}, {1, 0, 0}, {1, 1, 1}}};

/**
 * This class generates points uniformly over a tetrahedron in the form of
 * four-component barycentric coordinates useable for arbitrary tetrahedra.
 * This generator works by obtaining uniformly distributed points
 * over the unit cube.
 * Then the points are identified as being inside one of five cube-filling
 * tetrahedra*:
 * @f[
 *   T_0 = ((1, 1, 1), (0, 0, 1), (1, 0, 0), (0, 1, 0))
 *   T_1 = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
 *   T_2 = ((1, 0, 1), (0, 0, 1), (1, 0, 0), (1, 1, 1))
 *   T_3 = ((0, 1, 1), (0, 1, 0), (0, 0, 1), (1, 1, 1))
 *   T_4 = ((1, 1, 0), (0, 1, 0), (1, 0, 0), (1, 1, 1))
 * @f]
 * Finally, each point is converted to the barycentric coordinates of its
 * containing tetrahedron and returned.
 *
 * This distribution uses the uniform real distribution over its default range
 * of @f$ [0, 1) @f$ and is thus stateless.
 *
 * This technique is a generalization of the procedure used in the uniform
 * inside triangle distribution.  Adapting this to higher dimensions is
 * expected to be very difficult and of diminishing utility.
 *
 * *This tetrahedralization can be constructed by clipping every other vertex
 * of the cube through the adjacent vertices.  This produces four tetrahedral
 * pyramids of base edge length @f$ \sqrt{2} @f$ and sloping edge length 1
 * and leaves behind a regular tetrahedron of edge length @f$ \sqrt{2} @f$.
 */
template<class RealTp = double>
  class uniform_inside_tetrahedron_distribution
  {
  public:

    using result_type = std::array<RealTp, 4>;
    struct param_type
    { };

    uniform_inside_tetrahedron_distribution()
    : _M_urd(std::numeric_limits<RealTp>::min(), RealTp{1})
    { }

    uniform_inside_tetrahedron_distribution(param_type)
    : _M_urd(std::numeric_limits<RealTp>::min(), RealTp{1})
    { }

    param_type
    param() const
    { return param_type(); }

    void
    param(param_type)
    { }

    void
    reset()
    { this->_M_urd.reset(); }

    template<typename Generator>
      result_type
      operator()(Generator& __gen)
      {
	while (true)
	{
	  //  Pick three random barycentric coordinates.
	  std::array<RealTp, 3> __point;
	  __point[0] = this->_M_urd(__gen);
	  __point[1] = this->_M_urd(__gen);
	  __point[2] = this->_M_urd(__gen);

          if (auto bary = barycenter(tetra0<RealTp>, __point); valid(bary))
	    return bary;
          if (auto bary = barycenter(tetra1<RealTp>, __point); valid(bary))
	    return bary;
          if (auto bary = barycenter(tetra2<RealTp>, __point); valid(bary))
	    return bary;
          if (auto bary = barycenter(tetra3<RealTp>, __point); valid(bary))
	    return bary;
          if (auto bary = barycenter(tetra4<RealTp>, __point); valid(bary))
	    return bary;
	}
      }

    template<typename Generator>
      result_type
      operator()(Generator& __gen, param_type)
      { return this->operator()(__gen); }

    constexpr result_type
    min() const
    { return std::array<RealTp, 3>{0, 0, 0, 0}; }

    constexpr result_type
    max() const
    { return std::array<RealTp, 3>{1, 1, 1, 1}; }

  private:

    std::uniform_real_distribution<RealTp> _M_urd;
  };

// This thing is really stateless.
template<typename RealTp>
  constexpr bool
  operator==(const uniform_inside_tetrahedron_distribution<RealTp>&,
	     const uniform_inside_tetrahedron_distribution<RealTp>&)
  { return true; }

template<typename RealTp>
  constexpr bool
  operator!=(const uniform_inside_tetrahedron_distribution<RealTp>&,
	     const uniform_inside_tetrahedron_distribution<RealTp>&)
  { return false; }

template<typename CharT, typename Traits, typename RealTp>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const uniform_inside_tetrahedron_distribution<RealTp>&)
  { return __os; }

template<typename CharT, typename Traits, typename RealTp>
  std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT, Traits>& __is,
	     const uniform_inside_tetrahedron_distribution<RealTp>&)
  { return __is; }

int
main()
{
  using array_t = uniform_inside_tetrahedron_distribution<>::result_type;
  using real_t = array_t::value_type;

  std::cout.precision(std::numeric_limits<real_t>::digits10);
  auto w = 6 + std::cout.precision();

  uniform_inside_tetrahedron_distribution<> uit;

  std::cout << "uit = " << uit << '\n';

  tetra<real_t> tetra0;
  tetra<real_t> tetra1({6, 6, 6}, {-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

  std::ofstream iv("tetras.iv");
  iv << "#Inventor V2.1 ascii\n";
  iv << "Separator {\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";

  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 0; i < 100000; ++i)
    {
      auto bary = uit(gen);
      //std::cout << bary << '\n';
      std::cout << ' ' << std::setw(w) << bary[0]
		<< ' ' << std::setw(w) << bary[1]
		<< ' ' << std::setw(w) << bary[2]
		<< ' ' << std::setw(w) << bary[3] << '\n';
      iv << tetra0(bary) << ",\n";
      iv << tetra1(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";
  iv << "}\n";
}

