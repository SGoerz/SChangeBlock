# require(robcp)
# 
# modifChol2 <- function(x, ...)
# {
#   res <- modifChol(x, ...)
#   swaps <- attr(res, "swaps") + 1
#   n <- nrow(x)
#   i <- 1:n
#   sapply(n:1, function(j) {temp <- i[swaps[j]]; i[swaps[j]] <<- i[j]; i[j] <<- temp})
#   return(res[, i])
# }

#' @export
bandwidth <- function(X, p1, p2)
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
}


#' @importFrom robcp modifChol
#' @export
decorr <- function(X, lags, method = 1, separable = FALSE)
{
  ## because of autocovariance definition resp. optional kernel function
  lags <- lags + 1
  x <- as.vector(X)

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
  
  # if(method == 1)
  # {
  #   res <- svd(M)
  #   M1sq <- res$u %*% diag(1 / sqrt(res$d)) %*% t(res$v)
  # } else if(method == 2)
  # {
  #   Mchol <- tryCatch(chol(M), error = function(e) e)
  #   if("error" %in% class(Mchol)) Mchol <- modifChol(M)
  #   
  #   M1sq <- solve(Mchol)
  # } else if(method == 3 || method == "ar")
  # {
  #   Y <- arima(x, order = lags[1:3], 
  #              seasonal = list(order = lags[4:6], period = nrow(X)))$residuals
  #   dim(Y) <- dim(X)
  #   attributes(Y)$tsp <- NULL
  #   Y <- Y * 100
  #   class(Y) <- "RandomField"
  #   return(Y)
  # } 
  
  y <- t(M.invsqrt) %*% (x - mean(x))
  Y <- matrix(y, ncol = ncol(X))
  class(Y) <- "RandomField"
  return(Y)
}

invsqrt <- function(M, method = "")
{
  if(method == 1)
  {
    Mchol <- tryCatch(chol(M), error = function(e) e)
    if("error" %in% class(Mchol)) Mchol <- modifChol(M)
    return(solve(Mchol))
  } else if(method == 2)
  {
    res <- svd(M)
    return(res$u %*% diag(1 / sqrt(res$d)) %*% t(res$v))
  } else 
  {
    stop("Wrong method")
  }
}
