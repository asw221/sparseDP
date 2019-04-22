

## Not exported for now



hs.kappa.cm <- function(y, sigma) {
  K <- function(y, sigma) {
    2 * sqrt(2) * sigma * gsl::dawson(y / (sigma * sqrt(2))) / y
  }
  2 * sigma * (sqrt(2) * (y^2 + sigma^2) *
               gsl::dawson(y / (sigma * sqrt(2))) - sigma * y) /
    (y^3 * K(y, sigma))
}




