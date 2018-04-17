#ifndef DIRICHLET_DISTRIBUTION_H
#define DIRICHLET_DISTRIBUTION_H 1

#pragma GCC system_header

namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{
//_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * @brief A Dirichlet distribution for random numbers.
   *
   * The formula for the Dirichlet probability density function is
   * @f[
   *     p(\overrightarrow{x}|\overrightarrow{\alpha}) =
   *       \frac{1}{B(\overrightarrow{\alpha})}
   *       \prod_{i=1}^{k}x_i^{\alpha_i-1}
   * @f]
   * where @f$ \overrightarrow{x}, 0 < x_i < 1, \sum_{i=1}^{k}x_i = 1 @f$
   * and @f$ \overrightarrow{\alpha}, \alpha_i > 0 @f$ are
   * vectors of dimension @f$ k @f$
   * and where
   * @f[
   *    B(\overrightarrow{\alpha})
   *       = \frac{\prod_{i=1}^{k}\Gamma(\alpha_i)}
   *              {\Gamma(\sum_{i=1}^{k}\alpha_i)}
   * @f]
   * is the multi-dimensional beta function.
   *
   * The Dirichlet distribution is a multi-variate beta distribution over
   * k-dimensional simplices.
   */
  template<std::size_t _Dim, typename _RealTp = double>
    class dirichlet_distribution
    {
      static_assert(std::is_floating_point<_RealTp>::value,
		    "template argument not a floating point type");
      static_assert(_Dim != 0, "dimension is zero");

    public:
      /** The type of the range of the distribution. */
      typedef std::array<_RealTp, _Dim> result_type;
      /** Parameter type. */
      class param_type
      {

      public:
	typedef dirichlet_distribution<_Dim, _RealTp> distribution_type;
	friend class dirichlet_distribution<_Dim, _RealTp>;

	param_type()
	{ std::fill(_M_alpha.begin(), _M_alpha.end(), _RealTp{1}); }

	template<typename _ForwardIterator1>
	  param_type(_ForwardIterator1 __alphabegin,
		     _ForwardIterator1 __alphaend)
	{
	  __glibcxx_function_requires(_ForwardIteratorConcept<
				      _ForwardIterator1>)
	  _GLIBCXX_DEBUG_ASSERT(std::distance(__alphabegin, __alphaend)
				<= _Dim);

	  _M_init(__alphabegin, __alphaend);
	}

	param_type(std::initializer_list<_RealTp> __alpha)
	{
	  _GLIBCXX_DEBUG_ASSERT(__alpha.size() <= _Dim);
	  _M_init(__alpha.begin(), __alpha.end());
	}

	std::array<_RealTp, _Dim>
	alpha() const
	{ return _M_alpha; }

	friend bool
	operator==(const param_type& __p1, const param_type& __p2)
	{ return __p1._M_alpha == __p2._M_alpha; }

	friend bool
	operator!=(const param_type& __p1, const param_type& __p2)
	{ return !(__p1 == __p2); }

      private:
	template <typename _InputIterator1>
	  void _M_init(_InputIterator1 __alphabegin,
		       _InputIterator1 __alphaend);

	std::array<_RealTp, _Dim> _M_alpha;
      };

    public:
      dirichlet_distribution()
      : _M_param(), _M_gd()
      { }

      template<typename _ForwardIterator1>
	dirichlet_distribution(_ForwardIterator1 __alphabegin,
			       _ForwardIterator1 __alphaend)
	: _M_param(__alphabegin, __alphaend),
	  _M_gd()
	{ }

      dirichlet_distribution(std::initializer_list<_RealTp> __alpha)
      : _M_param(__alpha), _M_gd()
      { }

      explicit
      dirichlet_distribution(const param_type& __p)
      : _M_param(__p), _M_gd()
      { }

      /**
       * @brief Resets the distribution state.
       */
      void
      reset()
      { _M_gd.reset(); }

      /**
       * @brief Returns the alpha parameters of the distribution.
       */
      result_type
      alpha() const
      { return _M_param.alpha(); }

      /**
       * @brief Returns the parameter set of the distribution.
       */
      param_type
      param() const
      { return _M_param; }

      /**
       * @brief Sets the parameter set of the distribution.
       * @param __param The new parameter set of the distribution.
       */
      void
      param(const param_type& __param)
      { _M_param = __param; }

      /**
       * @brief Returns the greatest lower bound value of the distribution.
       */
      result_type
      min() const
      { result_type __res;
	__res.fill(std::numeric_limits<_RealTp>::lowest());
	return __res; }

      /**
       * @brief Returns the least upper bound value of the distribution.
       */
      result_type
      max() const
      { result_type __res;
	__res.fill(std::numeric_limits<_RealTp>::max());
	return __res; }

      /**
       * @brief Generating functions.
       */
      template<typename _UniformRandomNumberGenerator>
	result_type
	operator()(_UniformRandomNumberGenerator& __urng)
	{ return this->operator()(__urng, _M_param); }

      template<typename _UniformRandomNumberGenerator>
	result_type
	operator()(_UniformRandomNumberGenerator& __urng,
		   const param_type& __p);

      template<typename _ForwardIterator,
	       typename _UniformRandomNumberGenerator>
	void
	__generate(_ForwardIterator __f, _ForwardIterator __t,
		   _UniformRandomNumberGenerator& __urng)
	{ return this->__generate_impl(__f, __t, __urng, _M_param); }

      template<typename _ForwardIterator,
	       typename _UniformRandomNumberGenerator>
	void
	__generate(_ForwardIterator __f, _ForwardIterator __t,
		   _UniformRandomNumberGenerator& __urng,
		   const param_type& __p)
	{ return this->__generate_impl(__f, __t, __urng, __p); }

      /**
       * @brief Return true if two multi-variant normal distributions have
       *        the same parameters and the sequences that would
       *        be generated are equal.
       */
      template<size_t _Dim1, typename _RealTp1>
	friend bool
	operator==(const
		   __gnu_cxx::dirichlet_distribution<_Dim1, _RealTp1>& __d1,
		   const
		   __gnu_cxx::dirichlet_distribution<_Dim1, _RealTp1>& __d2);

      /**
       * @brief Inserts a %dirichlet_distribution random number distribution
       * @p __x into the output stream @p __os.
       *
       * @param __os An output stream.
       * @param __x  A %dirichlet_distribution random number distribution.
       *
       * @returns The output stream with the state of @p __x inserted or in
       * an error state.
       */
      template<size_t _Dim1, typename _RealTp1,
	       typename _CharT, typename _Traits>
	friend std::basic_ostream<_CharT, _Traits>&
	operator<<(std::basic_ostream<_CharT, _Traits>& __os,
		   const
		   __gnu_cxx::dirichlet_distribution<_Dim1, _RealTp1>& __x);

      /**
       * @brief Extracts a %dirichlet_distribution random number distribution
       * @p __x from the input stream @p __is.
       *
       * @param __is An input stream.
       * @param __x  A %dirichlet_distribution random number generator engine.
       *
       * @returns The input stream with @p __x extracted or in an error
       *          state.
       */
      template<size_t _Dim1, typename _RealTp1,
	       typename _CharT, typename _Traits>
	friend std::basic_istream<_CharT, _Traits>&
	operator>>(std::basic_istream<_CharT, _Traits>& __is,
		   __gnu_cxx::dirichlet_distribution<_Dim1, _RealTp1>&
		   __x);

    private:
      template<typename _ForwardIterator,
	       typename _UniformRandomNumberGenerator>
	void
	__generate_impl(_ForwardIterator __f, _ForwardIterator __t,
			_UniformRandomNumberGenerator& __urng,
			const param_type& __p);

      param_type _M_param;
      std::gamma_distribution<_RealTp> _M_gd;
  };

  /**
   * @brief Return true if two multi-variate normal distributions are
   * different.
   */
  template<size_t _Dim, typename _RealTp>
    inline bool
    operator!=(const __gnu_cxx::dirichlet_distribution<_Dim, _RealTp>& __d1,
	       const __gnu_cxx::dirichlet_distribution<_Dim, _RealTp>& __d2)
    { return !(__d1 == __d2); }

//_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include "dirichlet_distribution.tcc"

#endif // DIRICHLET_DISTRIBUTION_H
