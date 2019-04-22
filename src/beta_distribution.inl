
#include <iostream>
#include <random>


// Constructors & Destructors --------------------------------------------------

template< class RealType >
beta_distribution<RealType>::param_type::param_type(RealType a, RealType b) {
  if (a <= 0 || b <= 0)
    throw std::logic_error("Beta distribution parameters must be > 0");
  _alpha = a;
  _beta = b;
};


template< class RealType >
beta_distribution<RealType>::beta_distribution(RealType a, RealType b) {
  param(param_type(a, b));
};


template< class RealType >
beta_distribution<RealType>::beta_distribution(const param_type &par) {
  param(par);
};




// Utility Functions/Operators -------------------------------------------------

template< class RealType >
template< class Generator >
RealType beta_distribution<RealType>::operator() (Generator &g) {
  double x, y;
  x = _Gamma_A_(g);
  y = _Gamma_B_(g);
  return (RealType)(x / (x + y));
};



// Friend Functions ------------------------------------------------------------

template< class RealType >
bool operator==(
  const typename beta_distribution<RealType>::param_type &lhs,
  const typename beta_distribution<RealType>::param_type &rhs
) {
  const bool alphaSame = lhs._alpha == rhs._alpha;
  const bool betaSame = lhs._beta == rhs._beta;
  return alphaSame && betaSame;
};


template< class RealType >
bool operator!=(
  const typename beta_distribution<RealType>::param_type &lhs,
  const typename beta_distribution<RealType>::param_type &rhs
) {
  return !(lhs == rhs);
};


template< class RealType >
bool operator==(
  const beta_distribution<RealType> &lhs,
  const beta_distribution<RealType> &rhs
) {
  const bool paramsSame = lhs._par == rhs._par;
  const bool internalsSame = (lhs._Gamma_A_ == rhs._Gamma_A_) &&
    (lhs._Gamma_B_ == rhs._Gamma_B_);
  return paramsSame && internalsSame;
};


template< class RealType >
bool operator!=(
  const beta_distribution<RealType> &lhs,
  const beta_distribution<RealType> &rhs
) {
  return !(lhs == rhs);
};



// template< class RealType >
template< class RealType, class CharT, class Traits >
std::basic_ostream<CharT, Traits>& operator<<(
  std::basic_ostream<CharT, Traits> &ost,
  const typename beta_distribution<RealType>::param_type &pt
) {
  ost << pt._alpha << ", " << pt._beta;
  return ost;
};



// template< class RealType >
template< class RealType, class CharT, class Traits >
std::basic_ostream<CharT, Traits>& operator<<(
  std::basic_ostream<CharT, Traits> &ost,
  const beta_distribution<RealType> &d
) {
  ost << "Beta(" << d.param() << ")\n";
  return ost;
};



// Getters ---------------------------------------------------------------------

template< class RealType >
RealType beta_distribution<RealType>::alpha() const {
  return _par.alpha();
};


template< class RealType >
RealType beta_distribution<RealType>::beta() const {
  return _par.beta();
};


template< class RealType >
RealType beta_distribution<RealType>::max() const {
  return (RealType)1.0;
};


template< class RealType >
RealType beta_distribution<RealType>::min() const {
  return (RealType)0.0;
};


template< class RealType >
typename beta_distribution<RealType>::param_type
beta_distribution<RealType>::param() const {
  return _par;
};


template< class RealType >
RealType beta_distribution<RealType>::param_type::alpha() const {
  return _alpha;
};


template< class RealType >
RealType beta_distribution<RealType>::param_type::beta() const {
  return _beta;
};


// Setters ---------------------------------------------------------------------

template< class RealType >
void beta_distribution<RealType>::param(const param_type &par) {
  _par = par;
  _Gamma_A_.param(
    std::gamma_distribution<double>::param_type((double)par.alpha(), 1.0)
  );
  _Gamma_B_.param(
    std::gamma_distribution<double>::param_type((double)par.beta(), 1.0)
  );
};


template< class RealType >
void beta_distribution<RealType>::reset() {
  _Gamma_A_.reset();
  _Gamma_B_.reset();
};

