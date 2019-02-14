// /home/ed/bin/bin/g++ -std=gnu++17 -g -Iinclude -o test_von_mises_fisher_distribution test_von_mises_fisher_distribution.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <experimental/array>

#include <ext/von_mises_fisher_distribution.h>

template<typename Tp, size_t Dim>
  std::array<Tp, Dim>
  operator+(const std::array<Tp, Dim>& a, const std::array<Tp, Dim>& b)
  {
    std::array<Tp, Dim> c(a);
    for (unsigned int i = 0; i < Dim; ++i)
      c[i] += b[i];
    return c;
  }

template<typename Tp, size_t Dim>
  std::array<Tp, Dim>&
  operator+=(std::array<Tp, Dim>& a, const std::array<Tp, Dim>& b)
  {
    for (unsigned int i = 0; i < Dim; ++i)
      a[i] += b[i];
    return a;
  }

template<typename Tp, size_t Dim, typename Up>
  std::array<Tp, Dim>
  operator/(const std::array<Tp, Dim>& a, Up b)
  {
    std::array<Tp, Dim> c(a);
    for (unsigned int i = 0; i < Dim; ++i)
      c[i] /= b;
    return c;
  }

template<typename Tp, size_t Dim, typename Up>
  std::array<Tp, Dim>&
  operator/=(std::array<Tp, Dim>& a, Up b)
  {
    for (unsigned int i = 0; i < Dim; ++i)
      a[i] /= b;
    return a;
  }

template<typename Tp, size_t Dim>
  Tp
  norm(const std::array<Tp, Dim>& a)
  {
    Tp n{};
    for (unsigned int i = 0; i < Dim; ++i)
      n += a[i] * a[i];
    return std::sqrt(n);
  }

template<typename Tp, size_t Dim>
  std::array<Tp, Dim>
  normal(const std::array<Tp, Dim>& a)
  {
    std::array<Tp, Dim> c(a);
    auto n = norm(c);
    for (unsigned int i = 0; i < Dim; ++i)
      c[i] /= n;
    return c;
  }

template<typename Tp, size_t Dim>
  std::array<Tp, Dim>&
  normalize(std::array<Tp, Dim>& a)
  {
    auto n = norm(a);
    for (unsigned int i = 0; i < Dim; ++i)
      a[i] /= n;
    return a;
  }

/**
 * Output a std::array.
 */
template<std::size_t _Dim, typename _RealType,
	 typename _CharT, typename _Traits>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	     const std::array<_RealType, _Dim>& __x)
  {
    for (auto __k = 0U; __k < _Dim; ++__k)
      {
	if (__k > 0)
          __os << ',' << ' ';
	__os << std::setw(12) << __x[__k];
      }
    return __os;
  }

int
main()
{
  std::mt19937 re;

  std::cout << "\n\n  Dimension 2...\n\n";

  double phi0 = __gnu_cxx::__math_constants<double>::__pi_quarter;
  std::array<double, 2> mu2 = std::experimental::make_array(std::cos(phi0),
							    std::sin(phi0));
  std::array<decltype(mu2), 1> lambda2;
  __gnu_cxx::__detail::__make_basis(mu2, lambda2);
  std::cout << "  mu2       : " << mu2 << "  " << norm(mu2) << '\n';
  std::cout << "  lambda2[0]: " << lambda2[0]
            << "  " << norm(lambda2[0]) << '\n';

  __gnu_cxx::von_mises_fisher_distribution<2, double> vmd2(mu2, 100.0);
  const int num_samples2 = 1000;
  std::array<double, 2> mean2;
  for (auto i = 0; i < num_samples2; ++i)
    {
      auto dir2 = vmd2(re);
      mean2 += dir2;
      std::cout << " (" << dir2 << ") " << norm(dir2) << '\n';
    }
  mean2 /= num_samples2;
  std::cout << "  mean2    = " << mean2 << '\n';
  mean2 = normal(mean2);
  std::cout << "  meandir2 = " << mean2 << '\n';
  std::cout << "  mu2      = " << mu2 << '\n';
  std::cout << "  vmd2 = " << vmd2 << '\n';

  std::ofstream iv2("2d_rays.iv");
  iv2 << "#Inventor V2.1 ascii\n";
  iv2 << "Separator {\n";

  iv2 << "  BaseColor {\n";
  iv2 << "    rgb 0.0 0.0 1.0\n";
  iv2 << "  }\n";
  iv2 << "  Coordinate3 {\n";
  iv2 << "    point [\n";
  iv2 << "      0.0 0.0 0.0, " << mu2[0] << ' ' << mu2[1] << " 0.0\n";
  iv2 << "    ]\n";
  iv2 << "  }\n";
  iv2 << "  LineSet {}\n";

  iv2 << "  BaseColor {\n";
  iv2 << "    rgb 1.0 0.0 0.0\n";
  iv2 << "  }\n";
  iv2 << "  Coordinate3 {\n";
  iv2 << "    point [\n";
  iv2 << "      0.0 0.0 0.0, " << lambda2[0][0] << ' ' << lambda2[0][1] << " 0.0\n";
  iv2 << "    ]\n";
  iv2 << "  }\n";
  iv2 << "  LineSet {}\n";

  iv2 << "  BaseColor {\n";
  iv2 << "    rgb 1.0 1.0 0.0\n";
  iv2 << "  }\n";
  iv2 << "  Coordinate3 {\n";
  iv2 << "    point [\n";
  for (int i = 0; i < num_samples2; ++i)
    {
      auto dir2 = vmd2(re);
      iv2 << "      0.0 0.0 0.0, " << dir2[0] << ' ' << dir2[1] << " 0.0\n";
    }
  iv2 << "    ]\n";
  iv2 << "  }\n";
  iv2 << "  LineSet {}\n";

  iv2 << "}\n";


  std::cout << "\n\n  Dimension 3...\n\n";

  phi0 = __gnu_cxx::__math_constants<double>::__pi_quarter;
  double theta0 = __gnu_cxx::__math_constants<double>::__pi_third;
  std::array<double, 3> mu3 = std::experimental::make_array(std::cos(phi0) * std::sin(theta0),
							    std::sin(phi0) * std::sin(theta0),
							    std::cos(theta0));
  std::array<decltype(mu3), 2> lambda3;
  __gnu_cxx::__detail::__make_basis(mu3, lambda3);
  std::cout << "  mu3       : " << mu3 << "  " << norm(mu3) << "\n";
  for (auto k = 0; k < 2; ++k)
    std::cout << "  lambda3[" << k << "]: " << lambda3[k]
    	      << "  " << norm(lambda3[k]) << "\n";

  __gnu_cxx::von_mises_fisher_distribution<3, double> vmd3(mu3, 100.0);
  const int num_samples3 = 1000;
  std::array<double, 3> mean3;
  for (auto i = 0; i < num_samples3; ++i)
    {
      auto dir3 = vmd3(re);
      mean3 += dir3;
      std::cout << "  " << dir3 << "  " << norm(dir3) << '\n';
    }
  mean3 /= num_samples3;
  std::cout << "  mean3    = " << mean3 << '\n';
  mean3 = normal(mean3);
  std::cout << "  meandir3 = " << mean3 << '\n';
  std::cout << "  mu3      = " << mu3 << '\n';
  std::cout << "  vmd3 = " << vmd3 << '\n';

  std::ofstream iv3("3d_rays.iv");
  iv3 << "#Inventor V2.1 ascii\n";
  iv3 << "Separator {\n";

  iv3 << "  BaseColor {\n";
  iv3 << "    rgb 0.0 0.0 1.0\n";
  iv3 << "  }\n";
  iv3 << "  Coordinate3 {\n";
  iv3 << "    point [\n";
  iv3 << "      0.0 0.0 0.0, " << mu3[0] << ' ' << mu3[1] << ' ' << mu3[2] << '\n';
  iv3 << "    ]\n";
  iv3 << "  }\n";
  iv3 << "  LineSet {}\n";

  iv3 << "  BaseColor {\n";
  iv3 << "    rgb 1.0 0.0 0.0\n";
  iv3 << "  }\n";
  iv3 << "  Coordinate3 {\n";
  iv3 << "    point [\n";
  iv3 << "      0.0 0.0 0.0, " << lambda3[0][0] << ' ' << lambda3[0][1] << ' ' << lambda3[0][2] << '\n';
  iv3 << "    ]\n";
  iv3 << "  }\n";
  iv3 << "  LineSet {}\n";

  iv3 << "  BaseColor {\n";
  iv3 << "    rgb 0.0 1.0 0.0\n";
  iv3 << "  }\n";
  iv3 << "  Coordinate3 {\n";
  iv3 << "    point [\n";
  iv3 << "      0.0 0.0 0.0, " << lambda3[1][0] << ' ' << lambda3[1][1] << ' ' << lambda3[1][2] << '\n';
  iv3 << "    ]\n";
  iv3 << "  }\n";
  iv3 << "  LineSet {}\n";

  iv3 << "  BaseColor {\n";
  iv3 << "    rgb 1.0 1.0 0.0\n";
  iv3 << "  }\n";
  iv3 << "  Coordinate3 {\n";
  iv3 << "    point [\n";
  for (int i = 0; i < num_samples3; ++i)
    {
      auto dir3 = vmd3(re);
      iv3 << "      0.0 0.0 0.0, " << dir3[0] << ' ' << dir3[1] << ' ' << dir3[2] << '\n';
    }
  iv3 << "    ]\n";
  iv3 << "  }\n";
  iv3 << "  LineSet {}\n";

  iv3 << "}\n";


  std::cout << "\n\n  Dimension 4...\n\n";

  auto mu4 = std::experimental::make_array(0.1, 0.3, 0.4, std::sqrt(1.0 - 0.01 - 0.09 - 0.16));
  std::array<decltype(mu4), 3> lambda4;
  __gnu_cxx::__detail::__make_basis(mu4, lambda4);
  std::cout << "  mu4       : " << mu4 << "  " << norm(mu4) << '\n';
  for (auto k = 0; k < 3; ++k)
    std::cout << "  lambda4[" << k << "]: " << lambda4[k]
    	      << "  " << norm(lambda4[k]) << '\n';

  __gnu_cxx::von_mises_fisher_distribution<4, double> vmd4(mu4, 100.0);
  const int num_samples4 = 1000;
  std::array<double, 4> mean4;
  for (auto i = 0; i < num_samples4; ++i)
    {
      auto dir4 = vmd4(re);
      mean4 += dir4;
      std::cout << "  " << dir4 << "  " << norm(dir4) << '\n';
    }
  mean4 /= num_samples4;
  std::cout << "  mean4    = " << mean4 << '\n';
  mean4 = normal(mean4);
  std::cout << "  meandir4 = " << mean4 << '\n';
  std::cout << "  mu4      = " << mu4 << '\n';
  std::cout << "  vmd4 = " << vmd4 << '\n';
}
