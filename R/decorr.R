#'Bandwidth estimation 
#'
#'Calculate MSE-optimal bandwidths according to Andrews (1991).
#'
#'@param X numeric vector or matrix.
#'@param p1,p2 exponents for sample size n resp. estimated dependency, between 0 and 1.
#'
#'@details Bandwidth \eqn{\boldsymbol{b}^{(n,m)} = (b_1^{(n)}, b_2^{(m)})} is estimated via
#'         \deqn{b_i^{(k)} = \min\left(k-1, \max\left(1, k^{0.3} \left(\frac{2\rho_i}{1 - \rho_i^2}\right)^{0.3}\right)\right),} 
#'         where \eqn{\rho_1} and \eqn{\rho_2} are the mean row- resp. column-wise Spearman autocorrelations to lag 1.
#'
#'@returns A numeric vector containing one or two elements, depending on if a vector or matrix is supplied.
#'         In case of a matrix: the first value is the bandwidth for the row-wise 
#'         and the second one for the column-wise estimation.
#'         
#'@references Andrews, D. W. (1991). “Heteroskedasticity and autocorrelation consistent covariance 
#'            matrix estimation”. In: Econometrica: Journal of the Econometric Society, pp. 817–858.
#'
#'@examples 
#'X1 <- genField(c(50, 50), Theta = genTheta(1, 0.4))
#'bandwidth(X1, 0.3, 0.3)
#'
#'Theta <- matrix(c(0.08, 0.1, 0.08, 0.8, 1, 0.8, 0.08, 0.1, 0.08), ncol = 3)
#'X2 <- genField(c(50, 50), Theta = Theta)
#'bandwidth(X2, 1/3, 2/3)
#'
#' @export
bandwidth <- function(X, p1, p2)
{
  if(is.matrix(X))
  {
    n <- nrow(X)
    m <- ncol(X)
    
    rho1 <- mean(sapply(1:n, function(i) cor(X[i, -m], X[i, -1], method = "spearman")))
    rho2 <- mean(sapply(1:m, function(i) cor(X[-n, i], X[-1, i], method = "spearman")))
    
    param1 <- min(max(round(m^(p1) * ((2 * rho1) / (1 - rho1^2))^(p2)), 1), m-1)
    param2 <- min(max(round(n^(p1) * ((2 * rho2) / (1 - rho2^2))^(p2)), 1), n-1)
    
    if(is.na(param1)) param1 <- 1
    if(is.na(param2)) param2 <- 1
    
    return(c(param1, param2))
  } else
  {
    n <- length(X)
    rho <- abs(cor(X[1:(n-1)], X[2:n], method = "spearman"))
    return(min(max(round(n^p1 * (rho / (1 - rho))^p2), 1), n-1))
  }
}


#' De-correlation
#' 
#' @param X Random Field, numeric matrix, or time series
#' @param lags numeric vector containing two integer values: the bandwidths for the row- reps. column-wise autocovariance estimation.
#'             (Up to which lag should the autocovariances be estimated?) Lags must be smaller than the dimensions of X.
#' @param method 1L: square root of the matrix via [robcp::modifChol()], inversion via [solve()] \cr
#'               2L: square root and inversion via singular value decomposition \cr
#'               3L: square root via [expm::sqrtm()], inversion via [solve()]
#' @param separable if the autocovariance function is (assumed to be) separable in the two directions of X, 
#'                  those two autocovariances can be estimated separately and then combined as a Kronecker product.
#'
#' @importFrom robcp modifChol
#' @importFrom expm sqrtm
#' @export
decorr <- function(X, lags, method = 1L, separable = FALSE)
{
  if(is.ts(X)) 
  {
    tspX <- tsp(X)
    clsX <- NULL
  } else
  { 
    tspX <- NULL
    clsX <- class(X)
  }
  
  if(length(lags) == 1)
  {
    if(is.vector(X) | is.null(dim(X)))
    {
      lags[2] <- 0
      X <- as.matrix(X)
    } else
    {
      lags[2] <- lags[1]
    }
  }
  
  x <- as.vector(X)
  X <- as.matrix(X)

  if(separable)
  {
    M1 <- autocov(X, c(1, lags[2]), direction = 1)
    M2 <- autocov(X, c(lags[1], 1), direction = 2)
    
    M.invsqrt <- kronecker(invsqrt(M1, method), invsqrt(M2, method))
  } else 
  {
    M <- autocov(X, lags)
    M.invsqrt <- invsqrt(M, method)
  }

  y <- t(M.invsqrt) %*% (x - mean(x))
  if(ncol(X) > 1) Y <- matrix(y, ncol = ncol(X)) else Y <- as.vector(y)
  
  if(!is.null(tspX)) 
  {
    Y <- ts(Y)
    tsp(Y) <- tspX
  } else
  {
    class(Y) <- clsX
  }
  
  return(Y)
}

invsqrt <- function(M, method = 1L)
{
  if(method == 1L)
  {
    Mchol <- tryCatch(chol(M), error = function(e) e)
    if("error" %in% class(Mchol))
    {
      Mchol <- modifChol(M)
    }
    return(solve(Mchol))
  } else if(method == 2L)
  {
    res <- svd(M)
    return(res$u %*% diag(1 / sqrt(res$d)) %*% t(res$v))
  } else if(method == 3L)
  {
    Mchol <- tryCatch(chol(M), error = function(e) e)
    if("error" %in% class(Mchol))
    {
      Mchol <- sqrtm(M)
    }
    return(solve(Mchol))
  }
  {
    stop("Wrong method")
  }
}
