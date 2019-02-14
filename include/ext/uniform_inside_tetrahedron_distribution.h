#ifndef UNIFORM_INSIDE_TETRAHEDON_DISTRIBUTION_H
#define UNIFORM_INSIDE_TETRAHEDON_DISTRIBUTION_H 1

#pragma GCC system_header

#include <ext/random>
#include <ext/simplex.h>

namespace __gnu_test //_GLIBCXX_VISIBILITY(default)
{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

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

            if (auto bary = tetra0<RealTp>.barycenter(__point);
		bary_in_simplex(bary))
	      return bary;
            if (auto bary = tetra1<RealTp>.barycenter(__point);
		bary_in_simplex(bary))
	      return bary;
            if (auto bary = tetra2<RealTp>.barycenter(__point);
		bary_in_simplex(bary))
	      return bary;
            if (auto bary = tetra3<RealTp>.barycenter(__point);
		bary_in_simplex(bary))
	      return bary;
            if (auto bary = tetra4<RealTp>.barycenter(__point);
		bary_in_simplex(bary))
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
    inline std::basic_ostream<CharT,Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const uniform_inside_tetrahedron_distribution<RealTp>&)
    { return __os; }

  template<typename CharT, typename Traits, typename RealTp>
    inline std::basic_istream<CharT,Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is,
	       const uniform_inside_tetrahedron_distribution<RealTp>&)
    { return __is; }

//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_test

#endif // UNIFORM_INSIDE_TETRAHEDON_DISTRIBUTION_H
