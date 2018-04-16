#pragma GCC system_header

namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<std::size_t _Dimen, typename _RealType>
    template<typename _InputIterator1, typename _InputIterator2>
      void
      dirichlet_distribution<_Dimen, _RealType>::param_type::
      _M_init(_InputIterator1 __alphabegin, _InputIterator1 __alphaend)
      {
	__glibcxx_function_requires(_InputIteratorConcept<_InputIterator1>)
	std::fill(std::copy(__alphabegin, __alphaend, _M_alpha.begin()),
		  _M_alpha.end(), _RealType(0));
      }

  template<std::size_t _Dimen, typename _RealType>
    template<typename _UniformRandomNumberGenerator>
      typename dirichlet_distribution<_Dimen, _RealType>::result_type
      dirichlet_distribution<_Dimen, _RealType>::
      operator()(_UniformRandomNumberGenerator& __urng,
		 const param_type& __param)
      {
	result_type __ret;

	_M_gd.__generate(__ret.begin(), __ret.end(), __urng);

	for (size_t __i = _Dimen; __i > 0; --__i)
	  {
	    _RealType __sum = _RealType(0);
	    for (size_t __j = __i; __j > 0; --__j)
	      __sum += __ret[__j - 1];
	    __ret[__i - 1] = __sum;
	  }

	return __ret;
      }

  template<std::size_t _Dimen, typename _RealType>
    template<typename _ForwardIterator, typename _UniformRandomNumberGenerator>
      void
      dirichlet_distribution<_Dimen, _RealType>::
      __generate_impl(_ForwardIterator __f, _ForwardIterator __t,
		      _UniformRandomNumberGenerator& __urng,
		      const param_type& __param)
      {
	__glibcxx_function_requires(_Mutable_ForwardIteratorConcept<
				    _ForwardIterator>)
	while (__f != __t)
	  *__f++ = this->operator()(__urng, __param);
      }

  template<size_t _Dimen, typename _RealType>
    bool
    operator==(const __gnu_cxx::dirichlet_distribution<_Dimen, _RealType>&
	       __d1,
	       const __gnu_cxx::dirichlet_distribution<_Dimen, _RealType>&
	       __d2)
    {
      return __d1._M_param == __d2._M_param && __d1._M_gd == __d2._M_gd;
    }

  template<size_t _Dimen, typename _RealType, typename _CharT, typename _Traits>
    std::basic_ostream<_CharT, _Traits>&
    operator<<(std::basic_ostream<_CharT, _Traits>& __os,
	       const __gnu_cxx::dirichlet_distribution<_Dimen, _RealType>& __x)
    {
      typedef std::basic_ostream<_CharT, _Traits>  __ostream_type;
      typedef typename __ostream_type::ios_base    __ios_base;

      const typename __ios_base::fmtflags __flags = __os.flags();
      const _CharT __fill = __os.fill();
      const std::streamsize __precision = __os.precision();
      const _CharT __space = __os.widen(' ');
      __os.flags(__ios_base::scientific | __ios_base::left);
      __os.fill(__space);
      __os.precision(std::numeric_limits<_RealType>::max_digits10);

      auto __alpha = __x._M_param.alpha();
      for (auto __it : __alpha)
	__os << __it << __space;

      __os << __x._M_gd;

      __os.flags(__flags);
      __os.fill(__fill);
      __os.precision(__precision);
      return __os;
    }

  template<size_t _Dimen, typename _RealType, typename _CharT, typename _Traits>
    std::basic_istream<_CharT, _Traits>&
    operator>>(std::basic_istream<_CharT, _Traits>& __is,
	       __gnu_cxx::dirichlet_distribution<_Dimen, _RealType>& __x)
    {
      typedef std::basic_istream<_CharT, _Traits>  __istream_type;
      typedef typename __istream_type::ios_base    __ios_base;

      const typename __ios_base::fmtflags __flags = __is.flags();
      __is.flags(__ios_base::dec | __ios_base::skipws);

      std::array<_RealType, _Dimen> __alpha;
      for (auto& __it : __alpha)
	__is >> __it;

      __is >> __x._M_gd;

      __x.param(typename dirichlet_distribution<_Dimen, _RealType>::
		param_type(__alpha.begin(), __alpha.end()));

      __is.flags(__flags);
      return __is;
    }

//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx
