# function to cut one block and take the mean over the remaining data
# cutMean <- function(x, ln, ex = c(2, 2))
# {
#   if(!is.matrix(x)) x <- as.matrix(x)# stop("x is not a matrix")
#   
#   n <- dim(x)
#   
#   mean(x[1:(ln[1] - ex[1]), 1:(ln[2] - ex[2])])
# }

##' Cut-Off
##' 
##' This function returns a vector of the cut block means for a given time series or random field X.
##' 
##' @param X data, numeric vector or matrix
##' @param ln l block length. Integer vector of length 1 or 2, depending on the number of dimensions of X, with strictly positive entries.
##' @param ex number of blocks to cut off per dimension. Integer vector of length 1 or 2.
##'
##' @details For each block of size \code{ln[1]} \eqn{\cdot} \code{ln[2]}, the last \code{ex[1]} rows and the last \code{ex[2]} columns are deleted. Then, the mean is taken over the remaining data points.
##'
##' @return If X is a matrix, the return value is a numeric matrix of dimension \code{floor(n[1] / l[1]) x floor(n[2] / l[2])}, where \code{n = dim(X)}.
##'         Else, if X is a vector, the return value is a numeric vector of length \code{floor(n / l)}, where \code{n = length(X)}.
##'
##' @examples
##' X <- genField(c(40, 80), H = 100, type = 2)
##' M <- cutoff(X, l = c(4, 4))
##' 
##' plot(X)
##' image(M)
##'
##' @export
cutoff <- function(X, ln, ex = c(1, 1))
{
  classX <- class(X)
  
  if(is.matrix(X))
  {
    n <- dim(X)
    
    indexX <- Vectorize(function(x) seq(x, n[1], ln[1]))
    indexY <- Vectorize(function(x) seq(x, n[2], ln[2]))
    
    indexX <- sort(unlist(indexX(1:(ln[1] - ex[1]))))
    indexY <- sort(unlist(indexY(1:(ln[2] - ex[2]))))
    
    Y <- matrix(Mu(X[indexX, indexY], l = ln - ex, includeEdge = TRUE),
                nrow = ceiling(n[1] / ln[1]))
  } else
  {
    n <- length(X)
    
    indexX <- Vectorize(function(x) seq(x, n[1], ln[1]))
    indexX <- sort(unlist(indexX(1:(ln[1] - ex[1]))))
    
    Y <- Mu(X[indexX], l = ln - ex, includeEdge = TRUE)
  }
  
  class(Y) <- classX
  if(classX == "ts") Y <- ts(Y, start = start(X), frequency = frequency(X) / ln)
  
  return(Y)
}
