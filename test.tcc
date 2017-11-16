

template < typename AscendingSequence, typename CumulativeDist >
  test_result
  kolmogorov_test(AscendingSequence const & seq,
		  CumulativeDist const & cd)
  {
    KUBUX_ENFORCE(is_nondescending(seq));

    _Real k_upper = upper_kolmogorov_measure(seq, cd);
    _Real k_lower = lower_kolmogorov_measure(seq, cd);
    
    int upper_verdict = outcome(k_upper, cdist_kolmogorov());
    int lower_verdict = outcome(k_lower, cdist_kolmogorov());

    std::ostringstream msg;

    if (upper_verdict > 0)
      msg << "upper Kolmogorov measure " << k_upper << " too high, ";
    else if (upper_verdict < 0)
      msg << "upper Kolmogorov measure " << k_upper << " too small, ";
    else
      msg << "upper Kolmogorov measure " << k_upper << " ok, ";

    if (lower_verdict > 0)
      msg << "lower Kolmogorov measure " << k_lower << " too high";
    else if (lower_verdict < 0)
      msg << "lower Kolmogorov measure " << k_lower << " too small";
    else
      msg << "lower Kolmogorov measure " << k_lower << " ok";

    return (std::make_pair
	     ((upper_verdict == 0) && (lower_verdict == 0),
	       msg.str()));
  }

  // Section: The bucket test
  // ========================
  /*
    We divide the range into subranges and check whether each
    quantile gets its expected share.

    This is generic for real and integral RNGs.
  */


template<typename _Real, typename _UInt, typename UniformRNG, typename AscendingSequence >
  test_result
  generic_bucket_test(UniformRNG & urng,
		      AscendingSequence const & subdivision,
		      _UInt num_draws,
		      _Real eps = epsilon)
  {
    KUBUX_ENFORCE(is_ascending(subdivision));

    typedef typename UniformRNG::result_type value_type;
    typedef std::vector<_Real>             count_sequence;

    // determine the expected number of hits
    count_sequence expected;
    expected.push_back(subdivision[0] - urng.min() + is_discrete(urng));
    for (count_sequence::size_type i = 1; i < subdivision.size(); ++i)
      {
	expected.push_back(subdivision[i] - subdivision[i-1]);
	assert(expected.back() > 0);
      }
    expected.push_back(urng.max() - subdivision.back());
    assert(expected.back() > 0);
    for (count_sequence::size_type i = 0; i <= subdivision.size(); ++i)
      expected[i] *= _Real(num_draws) / size(urng);

    // throw the dices
    count_sequence observed (expected.size(), 0);
    while (num_draws-- > 0)
      {
	value_type value = urng();
	count_sequence::size_type index =
	  std::lower_bound(subdivision.begin(), subdivision.end(), value)
	  - subdivision.begin();
	++observed[ index ];
      }

    // compare observed and expected
    _Real chi_squared = get_chi_squared(observed, expected);
    int verdict = outcome(chi_squared, cdist_chi_squared(expected.size()-1));

    std::ostringstream msg;
    msg << "chi-square " << chi_squared;
    if (verdict > 0)
      msg << " too high";
    else if (verdict < 0)
      msg << " too small";
    else
      msg << " ok";

    return (std::make_pair(verdict == 0, msg.str()));
  }
  
template<typename UniformRNG >
  test_result
  full_bucket_test(UniformRNG & urng,
		   _UInt num_draws,
		   _Real eps = epsilon)
  {
    /*
      This version of the bucket test is for integral
      RNGs with small range. We just have a bucket for
      each possible value.
    */
    typedef typename UniformRNG::result_type value_type;
    typedef std::vector<value_type >        value_sequence;
    value_sequence subdivision;
    for (value_type i = urng.min(); i < urng.max(); ++i)
      subdivision.push_back(i);
    return (generic_bucket_test(urng, subdivision, num_draws, eps));
  }
  
  
template<typename _UInt, typename UniformRNG >
  test_result
  random_bucket_test(UniformRNG & urng,
		     _UInt num_buckets,
		     _UInt num_draws,
		     _Real eps = epsilon)
  {
    /*
      This version of the bucket test chooses the
      specified number of buckets randomly (using the
      universal_rng).
    */
    KUBUX_ENFORCE(num_buckets > 1);
    KUBUX_ENFORCE(num_buckets <= size(urng) || ! is_discrete(urng));
    
    typedef typename UniformRNG::result_type value_type;
    typedef std::vector<value_type >        value_sequence;

    // create the subdivision
    value_sequence subdivision;
    subdivision.reserve(num_buckets - 1);
    universal_rng the_rng;
    while (subdivision.size() != num_buckets -1)
      {
	for (_UInt i =  subdivision.size() + 1;
	      i < num_buckets; ++i)
	  {
	    value_type candidate = the_rng(urng.min(), urng.max());
	    if (urng.min() < candidate && candidate < urng.max())
	      subdivision.push_back(candidate);
	  }
	std::sort(subdivision.begin(),
		   subdivision.end());
	typename value_sequence::iterator last =
	  std::unique(subdivision.begin(), subdivision.end());
	subdivision.erase(last, subdivision.end());
      }
    return (generic_bucket_test(urng, subdivision, num_draws, eps));
  }


template<typename _UInt, typename UniformRNG >
  test_result
  bucket_test(UniformRNG & urng,
	      _UInt num_draws,
	      _Real eps = epsilon)
  {
    if (is_discrete(urng) &&	size(urng) <= 30000)
      return (full_bucket_test(urng, num_draws, eps));
    else
      return (random_bucket_test(urng, 10000, num_draws, eps));
  }

  
  // Section: The congruence test
  // ============================
  /*
    We divide the range into congruence classes modulo
    a given modulus. Then we check whether each class gets
    the expected share.

    This test only applies to integral RNGs.

    This test is currently not used.
  */
  
template<typename _UInt, typename IntType >
  _UInt
  mod(_UInt modulus, IntType i)
  {
    return (((i % modulus) + modulus) % modulus);
  }
  
template<typename UniformRNG >
  test_result
  congruence_test(UniformRNG & urng,
		  _UInt num_draws,
		  _UInt modulus,
		  _Real eps = epsilon)
  {
    KUBUX_ENFORCE(modulus <= size(urng));
    KUBUX_ENFORCE(is_discrete(urng));
    
    // compile expected hit counts
    _Real range = size(urng);
    _Real class_size = floorl(range / modulus);
    _UInt num_excess_slots = range - class_size * modulus;
    assert(num_excess_slots < modulus);
    _Real hits_expected = (class_size / range) * num_draws;
    std::vector<_Real> expected (modulus, hits_expected);
    for (_UInt i = 0; i < num_excess_slots; ++i)
      expected[i] += _Real(num_draws) / range;
 
    // compile observed hit counts
    std::vector<_Real> counted (modulus, 0.0L);
    for (_UInt i = 0; i < num_draws; ++i)
      ++counted[ (urng() - urng.min()) % modulus ];
 
    // and evaluate
    _Real chi_squared = get_chi_squared(counted, expected);
    int verdict =
      outcome(chi_squared,
	       cdist_chi_squared(expected.size()-1));

    std::ostringstream msg;
    msg << "chi-square " << chi_squared;
    if (verdict > 0)
      msg << " too high";
    else if (verdict < 0)
      msg << " too small";
    else
      msg << " ok";
    
    return (std::make_pair(verdict == 0, msg.str()));
  }
				 

  /*
    Kolmogorov test for uniformity.
    This test directly applies only to the uniform
    distribution of a real random variable in [0,1).

    We also use it for discrete distributions of large
    range size simply by rescaling.
  */

template<typename _UInt, typename UniformRNG >
  test_result
  uniform_test(UniformRNG & urng,
	       _UInt num_draws,
	       _Real eps = epsilon)
  {
    _Real rng_size = size(urng);
    std::vector<_Real> real (num_draws, 0.0L);
    for (unsigned int i = 0; i < num_draws; ++i)
      {
	real[i] = _Real(urng()) - _Real(urng.min());
	real[i] /= rng_size;
      }

    std::sort(real.begin(), real.end());
    return (kolmogorov_test(real, cdist_normalized_uniform()));
  }


  /*
    Minimum gap size test.
    This test compares the empirical distribution of the
    minimum circula gap to the theoretical one.

    This test applies directly to uniform reals in [0,1).

    We also use it for discrete RNGs with very large range
    size (>=2^24) by simple rescaling.
  */

template<typename _UInt, typename UniformRNG >
  test_result
  minimum_gap_test(UniformRNG & urng,
		   _UInt num_draws,
		   _UInt num_runs = 500,
		   _Real eps = epsilon)
  {
    _Real rng_size = size(urng);
    
    std::vector<_Real> real;
    real.reserve(num_runs);
    while (real.size() < num_runs)
      {
	std::vector<_Real> draws;
	draws.reserve(num_draws);
	while (draws.size() < num_draws)
	  {
	    _Real draw = _Real(urng()) - _Real(urng.min());
	    draw /= rng_size;
	    assert(0.0L <= draw && draw <= 1.0L);
	    draws.push_back(draw);
	  }
	std::sort(draws.begin(), draws.end());
	_Real last_distance = 1.0L + draws.front() - draws.back();
	for (std::vector<_Real>::size_type i = 1; i < num_draws; ++i)
	  draws[i-1] = draws[i] - draws[i-1];
	draws.back() = last_distance;
	std::sort(draws.begin(), draws.end());
	real.push_back(draws.front());
      }

    std::sort(real.begin(), real.end());

    assert(real.size() == num_runs);
    assert(is_nondescending(real));
    assert(is_in_unit_interval(real));
    
    return (kolmogorov_test(real, cdist_minimum_gap(num_draws)));
  }
