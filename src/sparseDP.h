
#define _USE_MATH_DEFINES

#include <Eigen/Core>
#include <cmath>


#ifndef _SPARSE_DP_
#define _SPARSE_DP_

namespace sparseDP {

  using std::exp;
  using std::sqrt;

  
  template< class RealType = double >
  RealType gaussian_std_density(const RealType &x) {
    return exp(-0.5 * x * x) / sqrt(2 * M_PI);
  };

  template< class RealType = double >
  RealType gaussian_density(
    const RealType &x,
    const RealType &mean = 0.0,
    const RealType &sd = 1.0
  ) {
    return gaussian_std_density( (x - mean)/sd ) / sd;
  };


  
  template< class Distribution >
  void setTwoParameters(
    Distribution &d,
    typename Distribution::result_type a,
    typename Distribution::result_type b
  ) {
    d.param(typename Distribution::param_type(a, b));
  };
  

  // template< class RandomIt >
  // void scaleSumToOne( RandomIt first, RandomIt last ) {
  //   if (first >= last) {
  //     throw std::logic_error("scaleSumToOne : ptrs should be first, last");
  //   }
  //   else if (first + 1 == last) {
  //     *first /= *first;
  //   }
  //   else {
  //     auto sum = *first;
  //     for (RandomIt it = first + 1; it != last; it++)
  // 	sum += *it;
  //     for (RandomIt it = first; it != last; it++)
  // 	*it /= sum;
  //   }
  // };
  
};


#endif  // _SPARSE_DP_
