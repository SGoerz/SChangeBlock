#' Mixture distribution
#' 
#' Generates a random sample from a mixture normal distribution.
#'
#' @param n sample size, numeric.
#' @param q probability that observation is drawn from the contamination distribution, numeric.
#' @param h mean of the contamination distribution, numeric.
#' @param sigma standard deviation of the contamination distribution, numeric.
#'
#' @return Numeric vector of length n containing the random sample.
#'
#' @examples
#' # random sample with 0.01 chance of contamination distribution with mean 10
#' rmix(100)
#' 
#' # random sample with 0.01 chance of contamination distribution with standard deviation 10
#' # IMPORTANT: h needs to be set to 0!
#' rmix(100, h = 0, sigma = 1)
#' 
#' @export
rmix <- function(n, q = 0.01, h = 10, sigma = 1)
{
  x <- rnorm(n)
  index <- sample(n, rbinom(1, n, q))
  x[index] <- (x[index] + h) * sigma
  return(x)
}


#' ???????????????
#' 
#' @examples
#' blockOutliers(100)
#' 
#' @export
blockOutliers <- function(n, q = 0.1, h = 10)
{
  x <- rnorm(n)
  n1 <- floor(sqrt(n))
  n2 <- floor(sqrt(n1))
  # range <- sapply(1:n2, function(x) ((x - 1) * n1 + 1):((x - 1) * n1 + n2))
  range <- sapply(1:n2, function(x) ((x) * n1 - n2 + 1):((x) * n1))
  index <- sample(range, ceiling(n2^2 * q))
  x[index] <- x[index] + h
  x <- x - ceiling(n2^2 * q) * h / n
  return(x)
}

#' ???DELETE???
#'@export
distributions <- Vectorize(function(type = "rnorm")
{
  distr <- switch(type, rnorm = rnorm,
                  rt3 = function(n) rt(n, 3), rt6 = function(n) rt(n, 6),
                  rchisq2 = function(n) rchisq(n, 2), rchisq3 = function(n) rchisq(n, 3),
                  rmix = function(n) rmix(n), rmixVar = function(n) rmix(n, h = 0, sigma = 10),
                  blockOutliers = blockOutliers, 
                  rbeta = function(n) rbeta(n, 10, 0.1), 
                  arima = function(n) arima.sim(list(ar = c(0.5, 0.3)), n))
  attr(distr, "name") <- type

  return(distr)
})

#' 
#' locFunctions <- Vectorize(function(type = "mean")
#' {
#'   fun <- switch(type, mean = mean, median = median)
#'   attr(fun, "name") <- type
#'
#'   return(fun)
#' })

