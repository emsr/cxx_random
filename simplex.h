#ifndef SIMPLEX_H
#define SIMPLEX_H 1

#include <array>

template<typename CharT, typename Traits, typename RealType>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealType, 3>& __arr)
  { return __os << ' ' << __arr[0] << ' ' << __arr[1] << ' ' << __arr[2]; }

template<typename CharT, typename Traits, typename RealTp>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealTp, 4>& __arr)
  { return __os << ' ' << __arr[0] << ' ' << __arr[1] << ' ' << __arr[2] << ' ' << __arr[3]; }

template<typename RealTp>
  struct tri
  {
    std::array<std::array<RealTp, 3>, 3> vert;

    tri()
    : vert{{{1,0,0},{0,1,0},{0,0,1}}}
    { }

    tri(const std::array<RealTp, 3>& A,
	const std::array<RealTp, 3>& B,
	const std::array<RealTp, 3>& C)
    : vert{A, B, C}
    { }

    std::array<RealTp, 3>
    operator()(const std::array<RealTp, 3>& bary) const
    {
      return {
	bary[0] * vert[0][0] + bary[1] * vert[1][0] + bary[2] * vert[2][0],
	bary[0] * vert[0][1] + bary[1] * vert[1][1] + bary[2] * vert[2][1],
	bary[0] * vert[0][2] + bary[1] * vert[1][2] + bary[2] * vert[2][2]};
    }
  };

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

#endif // SIMPLEX_H
