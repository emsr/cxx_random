#ifndef UNIFORM_INSIDE_TRIANGLE_DISTRIBUTION_H
#define UNIFORM_INSIDE_TRIANGLE_DISTRIBUTION_H 1

#pragma GCC system_header

#include <ext/random>
#include <ext/simplex.h>

namespace __gnu_test //_GLIBCXX_VISIBILITY(default)
{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
    inline std::basic_ostream<CharT,Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const uniform_inside_triangle_distribution<_RealType>&)
    { return __os; }

  template<typename CharT, typename Traits, typename _RealType>
    inline std::basic_istream<CharT,Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is,
	       const uniform_inside_triangle_distribution<_RealType>&)
    { return __is; }

//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_test

#endif // UNIFORM_INSIDE_TRIANGLE_DISTRIBUTION_H
