
#include <iostream>
#include <random>

/*
Implements a Beta distribution with functionality modeled after the
XXX_distribution classes in the <random> library.

Example use:
--------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

#include "beta_distribution.h"

// clang++ -g -Wall -O0 -std=c++17 -o test_beta_distribution.exe test_beta_distribution.cpp

int main() {
  std::random_device rd{};
  std::mt19937 gen{rd()};
 
  beta_distribution<> d{0.5, 0.5};
  // If X ~ Beta(a, b), the expected value of X is a / (a + b
  // values near 0 or 1 are most likely here
 
  std::map<double, int> hist{};
  for (int n = 0; n < 10000; n++) {
    hist[std::round(d(gen) * 10) / 10.0]++;
  }
  for (auto p : hist) {
    printf ("%.1f ", p.first);
    std::cout << std::setw(2);
    std::cout << std::string(p.second / 200, '*') << "\n";
  }
}


Example output:
--------------------------------------------------------------------------------
0.0 *******
0.1 *****
0.2 ****
0.3 ***
0.4 ***
0.5 ***
0.6 ***
0.7 ***
0.8 ***
0.9 *****
1.0 *******
 */



#ifndef _BETA_DISTRIBUTION_
#define _BETA_DISTRIBUTION_

template< class RealType = double >
class beta_distribution {
  
public:
  typedef RealType result_type;

  class param_type {
  private:
    RealType _alpha;  // shape 1
    RealType _beta;   // shape 2
    // Beta distribution has first moment _alpha / (_alpha + _beta)

  public:
    typedef beta_distribution<RealType> distribution_type;

    explicit param_type(RealType a = 1.0, RealType b = 1.0);

    RealType alpha() const;
    RealType beta() const;

    friend bool operator== (const param_type &lhs, const param_type &rhs);
    friend bool operator!= (const param_type &lhs, const param_type &rhs);

    template< class CharT, class Traits >
    friend std::basic_ostream<CharT, Traits>& operator<< (
      std::basic_ostream<CharT, Traits> &ost,
      const param_type &pt
    );
  };



  explicit beta_distribution(RealType a = 1.0, RealType b = 1.0);
  explicit beta_distribution(const param_type &par);

  template< class Generator >
  RealType operator() (Generator &g);

  RealType alpha() const;
  RealType beta() const;
  RealType max() const;
  RealType min() const;
  param_type param() const;

  void param(const param_type &par);
  void reset();
  
  friend bool operator== (
    const beta_distribution<RealType> &lhs,
    const beta_distribution<RealType> &rhs
  );

  friend bool operator!= (
    const beta_distribution<RealType> &lhs,
    const beta_distribution<RealType> &rhs
  );

  template< class CharT, class Traits >
  friend std::basic_ostream<CharT, Traits>& operator<< (
    std::basic_ostream<CharT, Traits> &ost,
    const beta_distribution<RealType> &d
  );

  
private:
  param_type _par;

  std::gamma_distribution<double> _Gamma_A_;
  std::gamma_distribution<double> _Gamma_B_;
};



#include "beta_distribution.inl"

#endif  // _BETA_DISTRIBUTION_
