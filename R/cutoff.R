# function to cut one block and take the mean over the remaining data
cutMean <- function(x, ln)
{
  if(!is.matrix(x))x <- as.matrix(x)# stop("x is not a matrix")
  
  n <- dim(x)
  
  mean(x[-ln[1], -ln[2]])
}

##' Cut-Off
##' 
##' This function returns a vector of the cut block means for a given time series or random field X.
##' 
##' @param X data, numeric vector or matrix
##' @param ln l block length. Integer vector of length 1 or 2, depending on the number of dimensions of X, with strictly positive entries.
##'
##' @details For each block of size \code{ln[1]} \eqn{\cdot} \code{ln[2]}, the \code{ln[1]}-th row and the \code{ln[2]}-th column are deleted. Then, the mean is taken over the remaining data points.
##'
##' @return A numeric vector of length \code{floor(n[1] / l[1]) * floor(n[2] / l[2])}, \code{n = dim(X)}.
##'
##' @examples
##' X <- genField(c(40, 80), H = 100, type = 2)
##' M <- cutoff(X, l = c(4, 4))
##' 
##' plot(X)
##' image(M)
##'
##' @export
cutoff <- function(X, ln)
{
  classX <- class(X)
  Y <- mu(X, ln = ln, function(x) cutMean(x, ln))
  class(Y) <- classX
  return(Y)
}