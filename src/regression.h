
#include <Eigen/Core>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <random>
#include <vector>

#include "SparseBrokenStick.h"
#include "sparseDP.h"


#ifndef _REGRESSION_
#define _REGRESSION_


const double _EPS_ = 1e-8;



void sampleClusterLabels(
  std::mt19937 &urng,               // modified
  std::vector<int> &clusterLabels,  // modified
  SparseBrokenStick &_stick_
);

void sampleClusterLabels(
  std::mt19937 &urng,               // modified
  std::vector<int> &clusterLabels,  // modified
  const int &j,
  const std::vector<double> &kernelWeights
);



std::vector<double> gaussianKernelClusters(
  const int &j,
  const Eigen::MatrixXd &X,
  const Eigen::VectorXd &residuals,
  const double &sigmaSq,
  const std::vector<int> &clusterLabels,
  const SparseBrokenStick &stick
);




Rcpp::List linearSparseDPRegression(
  const Eigen::MatrixXd &X,
  // const Eigen::MatrixXd &Z,
  const Eigen::VectorXd &y,
  const double &tauScalePrior,
  const int &truncateStickBreaking,
  const int &nSave,
  const int &thin,
  const int &burnin,
  const int &seed
) {
  std::mt19937 _rng_(seed);
  SparseBrokenStick stick(_rng_, truncateStickBreaking, tauScalePrior);

  std::normal_distribution<double> _StdNormal_;
  std::gamma_distribution<double> _Gamma_;       // shape, scale parameters

  const int n = y.size();
  const int P = X.cols();
  // const int D = Z.cols();
  // const bool shouldUpdateAlpha = D != 0;
  Eigen::VectorXd beta, residuals;
  double alpha, deltaAlpha, sigmaSq, residSS;
  int mcmcIter = 0, saveCount = 0;

  // Initialize parameters
  std::vector<int> clusterLabels(P);
  sampleClusterLabels(_rng_, clusterLabels, stick);

  beta = stick.beta(clusterLabels);
  alpha = _StdNormal_(_rng_);
  sigmaSq = 1 / _Gamma_(_rng_);

  // Initialize storage
  Eigen::MatrixXd betaSamples(nSave, P);
  Eigen::VectorXd alphaSamples(nSave);
  Eigen::VectorXd sigmaSamples(nSave);
  Eigen::VectorXd tauSamples(nSave);

  // Rcpp::Rcout << "-- Iteration 0 ---------------\n"
  // 	      << "sigma^2 = " << sigmaSq << "\n"
  // 	      << "alpha = " << alpha << "\n";
  // stick.Rprint();
  // Rcpp::Rcout << "------------------------------\n\n";

  // Main Gibbs sampling loop
  while (saveCount < nSave) {
    residuals = (y - X * stick.beta(clusterLabels)) -
      Eigen::VectorXd::Constant(n, alpha);
    
    // Update stick breaking components, residuals
    stick.update(_rng_, residuals, X, clusterLabels, sigmaSq);

    // Update intercept alpha, residuals
    deltaAlpha = std::sqrt(sigmaSq / n) * _StdNormal_(_rng_) +
      residuals.sum() / n;
    alpha += deltaAlpha;
    residuals -= Eigen::VectorXd::Constant(n, deltaAlpha);

    // Update residual sigma^2
    residSS = residuals.transpose() * residuals;
    sparseDP::setTwoParameters(_Gamma_, 0.5 * (n + stick.size()),
      2.0 / (residSS + stick.sumAtomsSq() / stick.globalScale()));
    // sparseDP::setTwoParameters(_Gamma_, 0.5 + 0.5 * n, 2.0 / (residSS + 1));
    sigmaSq = 1 / _Gamma_(_rng_);
    
    // Update cluster labels
    for (int j = 0; j < P; j++) {
      sampleClusterLabels(
        _rng_, clusterLabels, j,
        gaussianKernelClusters(j, X, residuals, sigmaSq, clusterLabels, stick)
      );
    }
 
    
    // Save samples and increment counters
    if (mcmcIter >= burnin && mcmcIter % thin == 0) {
      alphaSamples(saveCount) = alpha;
      betaSamples.row(saveCount) = stick.beta(clusterLabels).transpose();
      sigmaSamples(saveCount) = std::sqrt(sigmaSq);
      tauSamples(saveCount) = std::sqrt(stick.globalScale());
      saveCount++;
    }
    mcmcIter++;
  }

  // include log posterior later
  return Rcpp::List::create(
	   Rcpp::Named("alpha") = alphaSamples,
	   Rcpp::Named("beta") = betaSamples,
	   Rcpp::Named("sigma") = sigmaSamples,
	   Rcpp::Named("tau") = tauSamples
	 );
};
				    








void sampleClusterLabels(
  std::mt19937 &urng,               // modified
  std::vector<int> &clusterLabels,  // modified
  SparseBrokenStick &_stick_
) {
  std::discrete_distribution<int>
    _Multinom_(_stick_.weights_begin(), _stick_.weights_end());
  for (std::vector<int>::iterator label = clusterLabels.begin();
       label != clusterLabels.end(); label++)
    *label = _Multinom_(urng);
};


void sampleClusterLabels(
  std::mt19937 &urng,               // modified
  std::vector<int> &clusterLabels,  // modified
  const int &j,
  const std::vector<double> &kernelWeights
) {
  std::discrete_distribution<int>
    _Multinom_(kernelWeights.cbegin(), kernelWeights.cend());
  clusterLabels[j] = _Multinom_(urng);
};



std::vector<double> gaussianKernelClusters(
  const int &j,
  const Eigen::MatrixXd &X,
  const Eigen::VectorXd &residuals,
  const double &sigmaSq,
  const std::vector<int> &clusterLabels,
  const SparseBrokenStick &stick
) {
  // No bounds checking on atoms/clusterLabels, so careful
  const double currentAtom = stick.atom(clusterLabels[j]);
  const Eigen::VectorXd yStar = residuals + currentAtom * X.col(j);
  const double precHat = X.col(j).transpose() * X.col(j) + _EPS_;
  const double betaHat = (X.col(j).transpose() * yStar)[0] / precHat;
  const double sdHat = std::sqrt(sigmaSq / precHat);
  std::vector<double> kernel(stick.size());
  for (int h = 0; h < stick.size(); h++)
    kernel[h] = sparseDP::gaussian_density(stick.atom(h), betaHat, sdHat) *
      stick.weight(h) + _EPS_;
  return kernel;
};


#endif  // _REGRESSION_







// -----------------------------------------------------------------------------

      // Rcpp::Rcout << "-- Iteration " << (mcmcIter + 1) << " ---------------\n"
      // 		  << "sigma^2 = " << sigmaSq << "\n"
      // 		  << "alpha = " << alpha << "\n"
      // 		  << "residSS = " << residSS << "\n";
      // stick.Rprint();
      // 
      // Rcpp::Rcout << "\nCluster Labels:\n";
      // for (int j = 0; j < P; j++) {
      // 	if (j == 0)
      // 	  Rcpp::Rcout << clusterLabels[j];
      // 	else
      // 	  Rcpp::Rcout << ", " << clusterLabels[j];
      // }
      // Rcpp::Rcout << "\n";
      // Rcpp::Rcout << "------------------------------\n\n";
