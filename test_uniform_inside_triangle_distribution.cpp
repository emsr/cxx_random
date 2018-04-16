/*
$HOME/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp
*/

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "simplex.h"

template<class _RealType = double>
  class uniform_inside_triangle_distribution
  {
  public:

    using result_type = std::array<_RealType, 3>;
    struct param_type
    { };

    uniform_inside_triangle_distribution()
    : _M_urd(std::numeric_limits<_RealType>::min(), _RealType{1})
    { }

    uniform_inside_triangle_distribution(param_type)
    : _M_urd(std::numeric_limits<_RealType>::min(), _RealType{1})
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
	//  Pick two random barycentric coordinates.
	result_type __bary;
	__bary[1] = this->_M_urd(__gen);
	__bary[2] = this->_M_urd(__gen);
	if (__bary[1] + __bary[2] > _RealType{1})
	  {
	    __bary[1] = _RealType{1} - __bary[1];
	    __bary[2] = _RealType{1} - __bary[2];
	  }
	__bary[0] = _RealType{1} - __bary[1] - __bary[2];
	return __bary;
      }

    template<typename Generator>
      result_type
      operator()(Generator& __gen, param_type)
      { return this->operator()(__gen); }

    constexpr result_type
    min() const
    { return std::array<_RealType, 3>{0, 0, 0}; }

    constexpr result_type
    max() const
    { return std::array<_RealType, 3>{1, 1, 1}; }

  private:

    std::uniform_real_distribution<_RealType> _M_urd;
  };

// This thing is really stateless.
template<typename _RealType>
  bool
  operator==(const uniform_inside_triangle_distribution<_RealType>&,
	     const uniform_inside_triangle_distribution<_RealType>&)
  { return true; }

template<typename _RealType>
  bool
  operator!=(const uniform_inside_triangle_distribution<_RealType>&,
	     const uniform_inside_triangle_distribution<_RealType>&)
  { return false; }

template<typename CharT, typename Traits, typename _RealType>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const uniform_inside_triangle_distribution<_RealType>&)
  { return __os; }

template<typename CharT, typename Traits, typename _RealType>
  std::basic_istream<CharT,Traits>&
  operator>>(std::basic_istream<CharT, Traits>& __is,
	     const uniform_inside_triangle_distribution<_RealType>&)
  { return __is; }



/* Test Code */

int
main()
{
  using array_t = uniform_inside_triangle_distribution<>::result_type;
  using real_t = array_t::value_type;

  std::cout.precision(std::numeric_limits<real_t>::digits10);
  auto w = 6 + std::cout.precision();

  uniform_inside_triangle_distribution<> uit;

  std::cout << "uit = " << uit << '\n';

  tri<real_t> tri0;
  tri<real_t> tri1({-1, 3, 5}, {7, 2, -2}, {3, 3, -3});

  std::ofstream iv("tris.iv");
  iv << "#Inventor V2.1 ascii\n";
  iv << "Separator {\n";
  iv << "  Coordinate3 {\n";
  iv << "    point [\n";

  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 0; i < 1000; ++i)
    {
      auto bary = uit(gen);
      //std::cout << bary << '\n';
      std::cout << ' ' << std::setw(w) << bary[0]
		<< ' ' << std::setw(w) << bary[1]
		<< ' ' << std::setw(w) << bary[2] << '\n';
      iv << tri0(bary) << ",\n";
      iv << tri1(bary) << ",\n";
    }
  iv << "    ]\n";
  iv << "  }\n";
  iv << "  PointSet {}\n";
  iv << "}\n";
}

