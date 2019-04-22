
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "beta_distribution.h"


#ifndef _SPARSE_BROKEN_STICK_
#define _SPARSE_BROKEN_STICK_

//  - Blocked Gibbs sampler
//  - Horseshoe prior on stick-breaking atoms, \theta

class SparseBrokenStick {
private:
  typedef beta_distribution<double>::param_type _betaParam;
  typedef std::gamma_distribution<double>::param_type _gammaParam;  // shape, scale
  // typedef Eigen::Triplet<double> _Trips;
  
  double _BSq;
  double _eps_;
  double _precision;
  double _stickAlpha;
  double _sumAtomsSq;
  
  std::vector<double> _atoms;
  std::vector<double> _lambdaSq;
  std::vector<double> _weights;
  // std::vector<int> _clustCounts;

  beta_distribution<double> _Beta_;
  std::gamma_distribution<double> _Gamma_;
  std::normal_distribution<double> _StdNormal_;

  // Eigen::Triplet<double> getClusterTriplets(
  //   const std::vector<int> &clusterLabels,
  //   const int &cluster
  // );

  void setGammaParams(const double &shape, const double &rate);
  
public:

  explicit SparseBrokenStick(
    std::mt19937 &urng,
    int N = 25,
    double priorVar = 1.0,
    double alpha = 1.0,
    double B = 1.0
  );

  void update(
    std::mt19937 &urng,
    Eigen::VectorXd &residuals,
    const Eigen::MatrixXd &X,
    const std::vector<int> &clusterLabels,
    const double &sigmaSq
  );
  
  void updateWeights(std::mt19937 &urng, const std::vector<int> &clusterLabels);

  Eigen::VectorXd beta(const std::vector<int> &clusterLabels);
  // Eigen::SparseVector betaCluster(
  //   const std::vector<int> &clusterLabels,
  //   const int &cluster
  // );
  // Eigen::SparseVector betaCluster(
  //   const std::vector<int> &clusterLabels,
  //   const std::vector<int> &clusters
  // );

  double& operator[] (const int &index);

  int size() const;
  double atom(const int &index) const;
  double weight(const int &index) const;
  double globalScale() const;
  double sumAtomsSq() const;
  
  std::vector<double>::const_iterator cbegin() const;
  std::vector<double>::const_iterator cend() const;

  std::vector<double>::iterator weights_begin();
  std::vector<double>::iterator weights_end();

  // friend std::ostream operator<< (
  //   std::ostream &ost,
  //   const SparseBrokenStick &sbs
  // );
  void Rprint();
  
};




// -----------------------------------------------------------------------------

SparseBrokenStick::SparseBrokenStick(
  std::mt19937& urng,
  int N,
  double priorVar,
  double alpha,
  double B
) {
  if (N <= 0)
    throw std::logic_error("Stick-breaking must be truncated at N > 0");
  if (priorVar <= 0)
    throw std::logic_error("Prior variance must be > 0");
  if (alpha <= 0)
    throw std::logic_error("Stick-breaking alpha must be > 0");

  const double priorSd = std::sqrt(priorVar);
  double invWeightsCumProd = 1.0;
  std::vector<double> vecInvWeightsCumProd(N);

  _eps_ = 1e-8;
  _stickAlpha = alpha;
  _precision = std::pow(_StdNormal_(urng), 2) * priorVar;
  _BSq = B * B;
  _Beta_.param(_betaParam(1.0, _stickAlpha));

  setGammaParams(1.0, 1.0 / _BSq);
  _atoms.push_back(priorSd * _StdNormal_(urng));
  _lambdaSq.push_back(_Gamma_(urng));
  _weights.push_back(_Beta_(urng));
  vecInvWeightsCumProd[0] = 1.0;
  
  for (int h = 1; h < N; h++) {
    _atoms.push_back(priorSd * _StdNormal_(urng));
    _lambdaSq.push_back(_Gamma_(urng));
    _weights.push_back(_Beta_(urng));
    invWeightsCumProd *= 1.0 - _weights[h - 1];
    vecInvWeightsCumProd[h] = invWeightsCumProd;
  }
  for (int h = 1; h < N; h++)
    _weights[h] *= vecInvWeightsCumProd[h];
};




void SparseBrokenStick::update(
  std::mt19937 &urng,
  Eigen::VectorXd &residuals,
  const Eigen::MatrixXd &X,
  const std::vector<int> &clusterLabels,
  const double &sigmaSq
) {
  const int n = residuals.size(), P = clusterLabels.size();
  bool clusterOccupied;
  double thetaHatH, precH, nuH;
  double xi, sumThetaSq = 0.0;
  Eigen::VectorXd residStar, XStar;
  
  for (int h = 0; h < _atoms.size(); h++) {
    // compute cluster specific residuals and X
    clusterOccupied = false;
    residStar = residuals;
    XStar = Eigen::VectorXd::Zero(n);
    for (int j = 0; j < P; j++) {
      if (clusterLabels[j] == h) {
	XStar += X.col(j);
	clusterOccupied = true;
      }
    }

    // update local _lambdaSq
    setGammaParams(1.0, 1 + 1 / _lambdaSq[h]);
    nuH = 1 / _Gamma_(urng);
    setGammaParams(1.0, 1 / nuH + _atoms[h] * _atoms[h] * _precision / (2 * sigmaSq));
    _lambdaSq[h] = 1 / _Gamma_(urng);
    // _lambdaSq[h] = 1.0;

    // update _atoms
    if (clusterOccupied) {
      residStar += _atoms[h] * XStar;
      
      precH = (XStar.transpose() * XStar)[0] + _precision / _lambdaSq[h] + _eps_;
      // precH = (XStar.transpose() * XStar)[0] + sigmaSq * _precision;
      thetaHatH = (XStar.transpose() * residStar)[0] / precH;
      _atoms[h] = std::sqrt(sigmaSq / precH) * _StdNormal_(urng) + thetaHatH;
      
      // update residuals
      residuals = residStar - _atoms[h] * XStar;
    }
    else {
      _atoms[h] = std::sqrt(sigmaSq * _lambdaSq[h] / _precision) * _StdNormal_(urng);
      // _atoms[h] = std::sqrt(1 / _precision) * _StdNormal_(urng);
    }

    // updae _sumAtomsSq
    sumThetaSq += _atoms[h] * _atoms[h] / _lambdaSq[h];
  }

  // Update global _precision
  _sumAtomsSq = sumThetaSq;
  setGammaParams(1.0, 1 / _BSq + _precision);
  xi = _Gamma_(urng);
  setGammaParams(0.5 * (size() + 1), 1 / xi + 0.5 * _sumAtomsSq / sigmaSq);
  // setGammaParams(0.5 * size() + 0.5, 0.5 * _sumAtomsSq + 0.5);
  _precision = _Gamma_(urng);

  // Update _weights
  updateWeights(urng, clusterLabels);
};





void SparseBrokenStick::updateWeights(
  std::mt19937 &urng,
  const std::vector<int> &clusterLabels
) {
  const int Nc = size();     // Stick breaking truncation
  const int P = clusterLabels.size();  // Total number of parameters
  std::vector<int> clusterCounts(Nc);
  clusterCounts.assign(Nc, 0);
  for (std::vector<int>::const_iterator it = clusterLabels.cbegin();
       it != clusterLabels.cend(); it++) {
    if (*it >= 0 && *it < clusterCounts.size())
      clusterCounts[*it]++;
  }

  double countCumSum = clusterCounts[0], invWeightsCumProd = 1.0;
  std::vector<double> vecInvWeightsCumProd(Nc);
  _betaParam _bp(1 + clusterCounts[0], _stickAlpha + P - countCumSum);
  _Beta_.param(_bp);
  
  _weights[0] = _Beta_(urng);
  vecInvWeightsCumProd[0] = invWeightsCumProd;
  
  for (int h = 1; h < Nc; h++) {
    countCumSum += clusterCounts[h];
    _bp = _betaParam(1 + clusterCounts[h], _stickAlpha + P - countCumSum);
    _Beta_.param(_bp);
    _weights[h] = _Beta_(urng);
    invWeightsCumProd *= 1.0 - _weights[h - 1];
    vecInvWeightsCumProd[h] = invWeightsCumProd;
  }
  for (int h = 1; h < Nc; h++)
    _weights[h] *= vecInvWeightsCumProd[h];
};





// -----------------------------------------------------------------------------






void SparseBrokenStick::setGammaParams(const double &shape, const double &rate) {
  _gammaParam gPar(shape, 1 / rate);
  _Gamma_.param(gPar);
};




// -----------------------------------------------------------------------------

double& SparseBrokenStick::operator[](const int &index) {
  return _atoms[index];
};



Eigen::VectorXd SparseBrokenStick::beta(const std::vector<int> &clusterLabels) {
  const int P = clusterLabels.size();
  Eigen::VectorXd coefVector(P);
  for (int j = 0; j < P; j++)
    coefVector[j] = _atoms[clusterLabels[j]];
  return coefVector;
};




int SparseBrokenStick::size() const {
  return _atoms.size();
};



double SparseBrokenStick::atom(const int &index) const {
  return _atoms[index];
};


double SparseBrokenStick::weight(const int &index) const {
  return _weights[index];
};


double SparseBrokenStick::globalScale() const {
  return 1 / _precision;
};


double SparseBrokenStick::sumAtomsSq() const {
  return _sumAtomsSq;
};



std::vector<double>::const_iterator SparseBrokenStick::cbegin() const {
  return _atoms.cend();
};


std::vector<double>::const_iterator SparseBrokenStick::cend() const {
  return _atoms.cend();
};



std::vector<double>::iterator SparseBrokenStick::weights_begin() {
  return _weights.begin();
};


std::vector<double>::iterator SparseBrokenStick::weights_end() {
  return _weights.end();
};






void SparseBrokenStick::Rprint() {
  const int P = size();
  Rcpp::Rcout << "Blocked sampler truncation: " << P << "\n"
	      << "tau = " << 1 / std::sqrt(_precision) << "\n"
	      << "Atoms: ";
  for (int h = 0; h < P; h++) {
    if (h == 0)
      Rcpp::Rcout << _atoms[h];
    else
      Rcpp::Rcout << ", " << _atoms[h];
  }
  Rcpp::Rcout << "\n";

  Rcpp::Rcout << "Weights: ";
  for (int h = 0; h < P; h++) {
    if (h == 0)
      Rcpp::Rcout << _weights[h];
    else
      Rcpp::Rcout << ", " << _weights[h];
  }
  Rcpp::Rcout << "\n";

  Rcpp::Rcout << "lambda: ";
  for (int h = 0; h < P; h++) {
    if (h == 0)
      Rcpp::Rcout << std::sqrt(_lambdaSq[h]);
    else
      Rcpp::Rcout << ", " << std::sqrt(_lambdaSq[h]);
  }
  Rcpp::Rcout << "\n";
};






// -----------------------------------------------------------------------------






// std::ostream operator<< (
//   std::ostream &ost,
//   const SparseBrokenStick &sbs
// ) {
//   const int P = sbs.size();
//   ost << "Blocked sampler truncation: " << P << "\n"
//       << "Atoms: ";
//   for (int h = 0; h < P; h++) {
//     if (h == 0)
//       ost << sbs._atoms[h];
//     else
//       ost << ", " << sbs._atoms[h];
//   }
//   ost << "\n";
//
//   ost << "Weights: ";
//   for (int h = 0; h < P; h++) {
//     if (h == 0)
//       ost << sbs._weights[h];
//     else
//       ost << ", " << sbs._weights[h];
//   }
//   ost << "\n";
//
//   ost << "lambda: ";
//   for (int h = 0; h < P; h++) {
//     if (h == 0)
//       ost << std::sqrt(sbs._lambdaSq[h]);
//     else
//       ost << ", " << std::sqrt(sbs._lambdaSq[h]);
//   }
//   ost << "\n";
//   return ost;
// };
//
// typename SparseBrokenStick::_Trips SparseBrokenStick::getClusterTriplets(
//   const std::vector<int> &clusterLabels,
//   const int &cluster
// ) {
//   if (cluster < 0 || cluster >= _atoms.size())
//     throw std::logic_error("Desired cluster outside truncation");
//   const int P = clusterLabels.size();
//   const double atomH = _atoms[cluster];
//   std::vector<_Trips> trips;
//   for (int j = 0; j < P; j++) {
//     if (clusterLabels[j] == cluster)
//       trips.push_back(_Trips(j, 1, atomH));
//   }
//   return trips;
// };



// Eigen::SparseVector SparseBrokenStick::betaCluster(
//   const std::vector<int> &clusterLabels,
//   const int &cluster
// ) {
//   _Trips trips = getClusterTriplets(clusterLabels, cluster);
//   Eigen::SparseVector coefVector(clusterLabels.size());
//   coefVector.setFromTriplets(trips);
//   return coefVector;
// };


// Eigen::SparseVector betaCluster(
//   const std::vector<int> &clusterLabels,
//   const std::vector<int> &clusters
// ) {
//   _Trips wholeTrips;
//   _Trips partialTrips;
//   Eigen::SparseVector coefVector(clusterLabels.size());
//   for (std::vector<int>::iterator it = clusters.begin();
//        it != clusters.end(); it++) {
//     partialTrips = getClusterTriplets(clusterLabels, *it);
//     wholeTrips.insert(wholeTrips.end(),
// 		      partialTrips.begin(), partialTrips.end());
//   }
//   coefVector.setFromTriplets(wholeTrips);
//   return coefVector;
// };




#endif  // _SPARSE_BROKEN_STICK_
