template<typename _Real>
  class cdist_piecewise_constant
  {
    std::vector< _Real > the_intervals;
    std::vector< _Real > the_slopes;
    std::vector< _Real > the_sums;

    typedef std::vector< _Real >::size_type size_type;

  public:

    template < typename InputIteratorI, typename InputIteratorW >
    cdist_piecewise_constant ( InputIteratorI fromI, InputIteratorI toI,
			       InputIteratorW fromW )
    : the_intervals ( fromI, toI ),
      the_slopes ( the_intervals.size() - 1, 0 ),
      the_sums ( the_slopes.size(), 0 )
    {
      _Real the_total = 0;
      for ( size_type i = 0; i < the_slopes.size(); ++i ) {
	the_slopes[i] = *fromW;
	the_total += the_slopes[i];
	++fromW;
      }
      for ( size_type i = 0; i < the_slopes.size(); ++i ) {
	the_slopes[i] /= the_total;
      }
      for ( size_type i = 1; i < the_sums.size(); ++i ) {
	the_sums[i] = the_sums[i-1] + the_slopes[i-1];
      }
      for ( size_type i = 0; i < the_slopes.size(); ++i ) {
	the_slopes[i] /= the_intervals[i+1] - the_intervals[i];
      }
    }

    _Real operator() ( _Real x ) const {
      size_type index = lower_bound( the_intervals.begin(),
				     the_intervals.end(),
				     x )
	- the_intervals.begin();
      --index;
      if ( index >= the_slopes.size() ) {
	return ( 1 );
      }
      return ( the_sums[index] +
	       the_slopes[index]
	       *
	       ( x - the_intervals[index] ) );
    }
  };

template<typename _Real>
  class cdist_normalized_uniform
  {
  public:

    _Real operator() ( _Real x ) const {
      if ( x < 0.0L ) {
	return ( 0.0L );
      }
      if ( x > 1.0L ) {
	return ( 1.0L );
      }
      return ( x );
    }

  };

// The normal distribution
template<typename _Real>
  class cdist_normal {

    _Real mean;
    _Real stdev;

  public:

    cdist_normal ( _Real m = 0.0L, _Real d = 1.0L )
      : mean ( m )
      , stdev ( d )
    {}

    _Real operator() ( _Real x ) const {
      _Real const sqrt2 = 1.41421356237309504880L;
      return
	( 0.5L * ( 1.0L + erf( ( x - mean ) / ( stdev * sqrt2 ) ) ) );
    }
      
  };

template<typename _Real>
  class density_normal
  {
    
    _Real mean;
    _Real stdev;

  public:

    density_normal ( _Real m = 0.0L, _Real d = 1.0L )
      : mean ( m )
      , stdev ( d )
    {}

    _Real operator() ( _Real x ) const {
      x -= mean;
      x /= stdev;
      x *= x;
      x /= 2.0L;
      return ( exp( -x ) / ( stdev * sqrt2pi ) );
    }
      
  };

// The lognormal distribution
template<typename _Real>
  class cdist_lognormal
  {

    cdist_normal norm;

  public:

    cdist_lognormal ( _Real m = 0.0L, _Real d = 1.0L )
      : norm( m, d )
    {}

    _Real operator() ( _Real x ) const {
      return ( norm( log(x) ) );
    }
      
  };

// The exponential distribution
template<typename _Real>
  class cdist_exponential
  {

    _Real lambda;

  public:

    cdist_exponential ( _Real l = 1.0L )
      : lambda ( l ) 
    {}

    _Real operator() ( _Real x ) const {
      if ( x < 0.0L ) {
	return ( 0 );
      }
      return ( 1.0L - exp( -x/lambda ) );
    }
    
  };

// The Weibull distribution
template<typename _Real>
  class cdist_weibull
  {

    _Real k;
    _Real l;

  public:

    cdist_weibull ( _Real k_par, _Real l_par )
      : k ( k_par )
      , l ( l_par )
    {}

    _Real operator() ( _Real x ) const {
      return ( 1.0 - std::exp( - pow( x / l, k  ) ) );
    }

  };

// The extreme value distribution
template<typename _Real>
  class cdist_extreme_value
  {

    _Real a;
    _Real b;

  public:

    cdist_extreme_value ( _Real a_par, _Real b_par )
      : a ( a_par )
      , b ( b_par )
    {}

    _Real operator() ( _Real x ) const {
      return ( std::exp( - std::exp( - ( x - a ) / b ) ) );
    }

  };

// The chi-square distribution
template<typename _Real>
  class cdist_chi_squared
  {
    _Real k;

  public:
    
    cdist_chi_squared ( unsigned int val )
      : k ( val )
    {}

    _Real operator() ( _Real x ) const {
      return ( regularized_lower_gamma( k/2.0L, x/2.0L ) );
    }

  };

template <typename _Real, typename ObservedSequence, typename ExpectedSequence >
  _Real
  get_chi_squared( ObservedSequence const & obs,
		  ExpectedSequence const & exp )
  {
    assert( obs.size() == exp.size() );
    _Real chi_squared = 0.0L;
    for ( typename ObservedSequence::size_type i = 0;
	  i < obs.size(); ++i )
      chi_squared += sqr( obs[i] - exp[i] ) / exp[i];
    return ( chi_squared );
  }

// The Cauchy distribution
template<typename _Real>
  class cdist_cauchy
  {

    _Real a;
    _Real b;

  public:

    cdist_cauchy ( _Real a_par, _Real b_par )
      : a ( a_par )
      , b ( b_par )
    {}

    _Real operator() ( _Real x ) const {
      return  ( atanl( (x-a) / b ) / pi + 0.5L );
    }

  };


// The Gamma distribution
template<typename _Real>
  class cdist_gamma
  {

    _Real alpha;
    _Real beta;

  public:

    cdist_gamma ( _Real a, _Real b )
      : alpha ( a )
      , beta ( b )
    {}

    _Real operator() ( _Real x ) const {
      return ( incomplete_lower_gamma( alpha, x/beta )
	       /
	       gamma( alpha ) );
    }
    
  };


// The beta distribtion
template<typename _Real>
  class cdist_beta
  {

    _Real p;
    _Real q;

  public:

    cdist_beta ( _Real p_par, _Real q_par )
      : p ( p_par )
      , q ( q_par )
    {}

    _Real operator() ( _Real x ) const {
      if ( x < 0.0L ) { return 0.0L; }
      if ( x > 1.0L ) { return 1.0L; }
      return ( regularized_beta( x, p, q ) );
    }
    
  };

  
// The Fisher F distribtion
template<typename _Real>
  class cdist_fisher_f
  {

    _Real m;
    _Real n;

  public:

    cdist_fisher_f ( _Real m_par, _Real n_par )
      : m ( m_par )
      , n ( n_par )
    {}

    _Real operator() ( _Real x ) const {
      if ( x < 0.0L ) { return 0.0L; }
      return ( regularized_beta( 1.0L - n/(m*x+n), m/2, n/2 ) );
    }
    
  };


// Student's t distribution
template<typename _Real>
  class cdist_student_t
  {

    _Real n;

  public:

    cdist_student_t ( uhuge n_par )
      : n ( n_par )
    {}

    _Real operator() ( _Real t ) const {
      _Real a = sqrt( t*t + n );
      _Real x = ( t + a ) / ( 2.0L * a );
      _Real n_half = n/2.0L;
      return ( regularized_beta( x, n_half, n_half ) );
    }
    
  };

  
// The Kolmogorov distribution
template<typename _Real>
  class cdist_kolmogorov
  {
  public:

    _Real operator() ( _Real x ) const {
      return ( 1.0L - exp( -2.0L * sqr(x) ) );
    }

  };

// 
template < typename _Real, typename AscendingSequence, typename CumulativeDist >
  _Real upper_kolmogorov_measure ( AscendingSequence const & seq,
				   CumulativeDist cd ) {
    KUBUX_ASSERT( is_nondescending( seq ) );
    _Real the_max = 0;
    _Real index = 1;
    for ( typename AscendingSequence::const_iterator iter
	    = seq.begin(); iter != seq.end(); ++iter ) {
      _Real term = index - cd( *iter ) * seq.size();
      the_max = std::max( the_max, term );
      ++ index;
    }
    return ( the_max / root( seq.size() ) );
  }

// 
template < typename _Real, typename AscendingSequence, typename CumulativeDist >
  _Real lower_kolmogorov_measure ( AscendingSequence const & seq,
				   CumulativeDist cd ) {
    KUBUX_ASSERT( is_nondescending( seq ) );
    _Real the_max = 0;
    _Real index = 0;
    for ( typename AscendingSequence::const_iterator iter
	    = seq.begin(); iter != seq.end(); ++iter )
      {
	_Real term = cd( *iter ) * seq.size() - index;
	the_max = std::max( the_max, term );
	++ index;
      }
    return ( the_max / root( seq.size() ) );
  }


  // The minimum gap distribution
  //
  //  This introduces the order statistic of gaps of
  //  n random points on a circle of unit perimeter.
template<typename _Real>
  class cdist_minimum_gap
  {

    _Real num_points;

  public:

    cdist_minimum_gap ( uhuge n )
      : num_points ( n )
    {}

    _Real operator() ( _Real x ) const {
      if ( x <= 0.0L )
	return ( 0 );
      else if ( num_points*x >= 1.0L )
	return ( 1.0L );
      else
        return ( 1.0L - std::pow( 1.0L - num_points*x,
				num_points-1 ) );
    }

  };
