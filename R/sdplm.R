
sdplm <- function(X, y, tau.scale.prior = 1.0,
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



