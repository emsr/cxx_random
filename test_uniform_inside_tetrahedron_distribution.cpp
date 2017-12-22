/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp
*/

#include <random>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>

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

template<class RealTp = double>
  RealTp
  volume(const std::array<std::array<RealTp, 3>, 4>& tetra)
  {
    std::array<RealTp, 3> u, v, w;
    for (auto c = 0; c < 3; ++c)
      {
	u[c] = tetra[1][c] - tetra[0][c];
	v[c] = tetra[2][c] - tetra[0][c];
	w[c] = tetra[3][c] - tetra[0][c];
      }
    return w[0] * (u[1] * v[2] - u[2] * v[1])
	 + w[1] * (u[2] * v[0] - u[0] * v[2])
	 + w[2] * (u[0] * v[1] - u[1] * v[0]);
  }

template<class RealTp = double>
  std::array<RealTp, 4>
  barycenter(const std::array<std::array<RealTp, 3>, 4>& tetra,
	     const std::array<RealTp, 3>& point)
  {
    auto vol = volume(tetra);
    std::array<RealTp, 4> bary;
    for (auto v = 0; v < 4; ++v)
      {
	auto tet = tetra;
	for (auto c = 0; c < 3; ++c)
	  tet[v][c] = point[c];
	bary[v] = volume(tet) / vol;
      }

    return bary;
  }

template<class RealTp = double>
  bool
  valid(const std::array<RealTp, 4>& bary)
  {
    const auto eps = RealTp{10} * std::numeric_limits<RealTp>::epsilon();
    auto sum = RealTp{0};
    for (int v = 0; v < 4; ++v)
      {
	if (bary[v] < RealTp{0} || bary[v] > RealTp{1})
	  return false;
        sum += bary[v];
      }

    return std::abs(sum - RealTp{1}) < eps;
  }

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
    : _M_urd()
    { }

    uniform_inside_tetrahedron_distribution(param_type)
    : _M_urd()
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



/* Test Code */

template<typename CharT, typename Traits, typename RealTp>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealTp, 3>& __arr)
  { return __os << ' ' << __arr[0] << ' ' << __arr[1] << ' ' << __arr[2]; }

template<typename CharT, typename Traits, typename RealTp>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealTp, 4>& __arr)
  { return __os << ' ' << __arr[0] << ' ' << __arr[1] << ' ' << __arr[2] << ' ' << __arr[3]; }

/**
 * A little tetrahedron class.
 */
template<typename RealTp>
  struct tetra
  {
    std::array<std::array<RealTp, 3>, 4> vert;

    /**
     * Default constructor.
     * Construct a tetrahedron with one vertex at the origin
     * and the other three at x = 1, y = 1, and z = 1 respectively.
     */
    tetra()
    : vert{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}}
    { }

    /**
     * Constructor taking four vertices.
     * Imagining A to be 'above' the triangle formed by BCD the vertices
     * BCD should traverse the triangle in a 'clockwise' direction to form
     * a positive-volume tetrahedron so that in the usual right hand rule
     * convention all the triangle normals face outward from the tetrahedron.
     */
    tetra(const std::array<RealTp, 3>& A,
	  const std::array<RealTp, 3>& B,
	  const std::array<RealTp, 3>& C,
	  const std::array<RealTp, 3>& D)
    : vert{A, B, C, D}
    { }

    /**
     * Given a point in barycentric coordinates return the point in R3.
     */
    std::array<RealTp, 3>
    operator()(const std::array<RealTp, 4>& bary) const
    {
      return {
	bary[0] * vert[0][0] + bary[1] * vert[1][0] + bary[2] * vert[2][0] + bary[3] * vert[3][0],
	bary[0] * vert[0][1] + bary[1] * vert[1][1] + bary[2] * vert[2][1] + bary[3] * vert[3][1],
	bary[0] * vert[0][2] + bary[1] * vert[1][2] + bary[2] * vert[2][2] + bary[3] * vert[3][2]};
    }

    /**
     * Given a point in R3 return the point in barycentric coordinates.
     */
    std::array<RealTp, 4>
    operator()(const std::array<RealTp, 3>& point)
    {
      auto vol = volume(vert);
      std::array<RealTp, 4> bary;
      for (auto v = 0; v < 4; ++v)
	{
	  auto tet = vert;
	  for (auto c = 0; c < 3; ++c)
	    tet[v][c] = point[c];
	  bary[v] = volume(tet) / vol;
	}

      return bary;
    }

    /**
     * Return the volume of the tetrahedron.
     */
    RealTp
    volume() const
    {
      std::array<RealTp, 3> u, v, w;
      for (auto c = 0; c < 3; ++c)
	{
	  u[c] = vert[1][c] - vert[0][c];
	  v[c] = vert[2][c] - vert[0][c];
	  w[c] = vert[3][c] - vert[0][c];
	}
      return w[0] * (u[1] * v[2] - u[2] * v[1])
	   + w[1] * (u[2] * v[0] - u[0] * v[2])
	   + w[2] * (u[0] * v[1] - u[1] * v[0]);
    }
  };

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

