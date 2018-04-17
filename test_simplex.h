
template<typename CharT, typename Traits, typename RealType>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealType, 3>& __arr)
  { return __os << ' ' << __arr[0] << ' ' << __arr[1] << ' ' << __arr[2]; }

template<typename CharT, typename Traits, typename RealTp>
  std::basic_ostream<CharT,Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& __os,
	     const std::array<RealTp, 4>& __arr)
  {
    return __os << ' ' << __arr[0] << ' ' << __arr[1]
		<< ' ' << __arr[2] << ' ' << __arr[3];
  }

/*

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
*/
