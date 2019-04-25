
## sdplm
## -----------------------------------------------------------------------------
#' @title Sparse Dirichlet process linear regression
#'
#' @description
#' Fits a linear model with a Horseshoe-centered Dirichlet
#' process prior on the regression coefficients using a blocked Gibbs
#' sampler. The prior is very flexible and results in approximately sparse
#' coefficient clustering from iteration to iteration. Each call to
#' \code{sdplm} results in one MCMC chain.
#'
#' @param X
#' Design matrix of predictors
#'
#' @param y
#' Vector of outcomes
#'
#' @param intercept
#' Whether or not to include a model intercept. This parameter is
#' intended to provide a similar interface to
#' \code{glmnet::\link[glmnet]{glmnet}},
#' but is not currently in used. My eventual intent is to allow
#' users to pass in a separate matrix of predictors, \code{Z},
#' with associated coefficients given uninformative priors rather than
#' the sparse Dirichlet process formulation
#'
#' @param tau.scale.prior
#' Hyperparameter for the global scale of the Dirichlet process atoms,
#' \eqn{\tau_0}. Default is 1
#'
#' @param sb.truncation
#' Truncation level for the stick breaking process/number of atoms
#' kept in memory. Default is 25
#'
#' @param nsave
#' MCMC parameter: number of posterior samples to return. Default is 1000
#'
#' @param thin
#' MCMC parameter: thinning factor so only 1 in \code{thin} posterior
#' samples are returned. Default is 5. The total number of MCMC iterations
#' is \code{burnin} + \code{nsave} * \code{thin}
#'
#' @param burnin
#' MCMC parameter: number of samples to discard from the start of the
#' chain
#'
#' @param .rng.seed
#' Random number generator seed for thread safety and consistency.
#' Default is set from the clock
#'
#' @return
#' An S3 object of class \code{sdplm} containing posterior samples
#' of each parameter
#'
###
## @examples, @seealso, @details, @references
##
sdplm <- function(X, y, intercept = TRUE, tau.scale.prior = 1.0,
                  sb.truncation = 25, nsave = 1000, thin = 5,
                  burnin = 5000, .rng.seed = as.integer(Sys.time())
                  ) {
  if (nrow(X) != length(y))
    stop ("Data dimension mismatch: nrow(X) != length(y)")
  if (nsave <= 0)
    stop ("MCMC parameters: nsave must be > 0")
  if (thin <= 0)
    stop ("MCMC parameters: chian thinning factor must be > 0")
  if (burnin < 0)
    stop ("MCMC parameters: burnin/warmup period must be >= 0")
  if (thin > 100)
    warning ("MCMC parameters: chain thinning factor is ", thin, " iterations.\n",
             "Chain may take a long time to run.")
  tau.scale.prior <- abs(tau.scale.prior)
  .rng.seed <- abs(as.integer(.rng.seed))

  structure(
    .Call("sdplm_Cpp", X, y, tau.scale.prior, sb.truncation, nsave,
          thin, burnin, .rng.seed, PACKAGE = "sparseDP"),
    class = "sdplm")
}





sigma.sdplm <- function(object, ...)  mean(object$sigma)

coef.sdplm <- function(object, ...) c(mean(object$alpha), colMeans(object$beta))



