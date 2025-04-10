#'Long run variance estimation bandwidth(s)
#'
#'Calculate optimal bandwidths for the long run variance (lrv) estimation of a 2-dim random field.
#'
#'@param X numeric matrix.
#'@param p1,p2 exponents for sample size n resp. estimated dependency, between 0 and 1.
#'
#'@details .................
#'
#'@returns A numeric vector containing two elements: the first value is the bandwidth for the row-wise 
#'         and the second one for the column-wise estimation.
#'
#'@examples 
#'X1 <- genField(c(50, 50), Phi = genPhi(1, 0.4))
#'lrvBandwidth(X1, 1/3, 2/3)
#'
#'Phi <- matrix(c(0.08, 0.1, 0.08, 0.8, 1, 0.8, 0.08, 0.1, 0.08), ncol = 3)
#'X2 <- genField(c(50, 50), Phi = Phi)
#'lrvBandwidth(X2, 1/3, 2/3)
#'
#'@export
lrvBandwidth <- function(X, p1, p2)
{
  n <- nrow(X)
  m <- ncol(X)
  k1 <- abs(mean(apply(X, 1, function(x) cor(x[1:(m-1)], x[2:m], method = "spearman"))))
  k2 <- abs(mean(apply(X, 2, function(x) cor(x[1:(n-1)], x[2:n], method = "spearman"))))
  
  b1 <- max(n^p1 * (k1 / (1 - k1))^p2, 1)
  b2 <- max(n^p1 * (k2 / (1 - k2))^p2, 1)
  
  return(c(b1 = b1, b2 = b2))
}

