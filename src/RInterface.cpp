
#include <Rcpp.h>
#include <RcppEigen.h>

#include "regression.h"


// --- Sparse Dirichlet Process Linear Model -----------------------------------
extern "C" SEXP sdplm_Cpp(
  const SEXP X_,
  const SEXP y_,
  const SEXP tauScalePrior_,
  const SEXP truncation_,
  const SEXP nSave_,
  const SEXP thin_,
  const SEXP burnin_,
  const SEXP seed_
) {
  try {
    const Eigen::Map<Eigen::MatrixXd>
      X(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X_));
    const Eigen::Map<Eigen::VectorXd>
      y(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(y_));
    const double tauScalePrior(Rcpp::as<double>(tauScalePrior_));
    const int truncation(Rcpp::as<int>(truncation_));
    const int nSave(Rcpp::as<int>(nSave_));
    const int thin(Rcpp::as<int>(thin_));
    const int burnin(Rcpp::as<int>(burnin_));
    const int seed(Rcpp::as<int>(seed_));

    return Rcpp::wrap(linearSparseDPRegression(
		        X, y, tauScalePrior, truncation, nSave,
			thin, burnin, seed)
		      );
  }
  catch (std::exception& _ex_) {
    forward_exception_to_r(_ex_);
  }
  catch (...) {
    ::Rf_error("C++ exception (unknown cause)");
  }
  return R_NilValue;  // not reached
};
  
