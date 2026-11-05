#' Build block index matrices
#' 
#' Builds a matrix (or vector) containing the indices which block the corresponding 
#' observation in a vector or matrix of the same dimension belongs to.
#' 
#' @param n sample size(s), integer (vector).
#' @param bn desired number of blocks (per dimension). Integer (vector).
#' 
#' @details
#' This function is especially useful if you cannot have equally sized blocks 
#' for performing the \code{block_test}. First it divides the observations in \code{bn}
#' equally sized blocks. Then it calculates the number of remaining observations
#' \eqn{b_n^{rest} = n - b_n \cdot \lfloor n / b_n \rfloor} and distributes them evenly 
#' among the first \eqn{b_n^{rest}} blocks. 
#' 
#' @returns Strictly positive integer vector or matrix having the dimensions \code{n}. 
#' 
#' @examples
#' X <- genField(c(11, 14))
#' group <- build_blocks(dim(X), c(2, 3))
#' group
#' 
#' Mu(X, group = group)
#' 
#' @export
build_blocks <- function(n, bn)
{
  if(length(n) == 1) n[2] <- 1
  if(length(bn) == 1) bn[2] <- 1
  
  ln <- floor(n / bn)

  bndiff <- n - ln * bn
  
  nblocks <- prod(bn)
  groups <- matrix(1:nblocks, nrow = bn[1])
  
  A <- kronecker(groups[0:bndiff[1], 0:bndiff[2], drop = FALSE],
                 matrix(1, nrow = ln[1] + 1, ncol = ln[2] + 1))
  B <- kronecker(matrix(groups[0:bndiff[1], (bndiff[2] + 1):bn[2]], nrow = bndiff[1]),
                 matrix(1, nrow = ln[1] + 1, ncol = ln[2]))
  C <- kronecker(groups[(bndiff[1] + 1):bn[1], 0:bndiff[2], drop = FALSE],
                 matrix(1, nrow = ln[1], ncol = ln[2] + 1))
  D <- kronecker(matrix(groups[(bndiff[1] + 1):bn[1], (bndiff[2] + 1):bn[2]], nrow = bn[1] - bndiff[1]),
                 matrix(1, nrow = ln[1], ncol = ln[2]))
  
  if(length(A) == 0) A <- NULL
  if(length(B) == 0) B <- NULL
  if(length(C) == 0) C <- NULL
  if(length(D) == 0) D <- NULL
  
  rbind(cbind(A, B), cbind(C, D))
}

#' Optimal parameter s
#' 
#' Calculates the best parameter \eqn{\tilde{s}} for a given approximation s, such that \eqn{n \; % \; \[n^{s}\] = 0}.
#' 
#' @param n sample size(s), integer (vector).
#' @param s desired exponent, \eqn{0.5 \leq s \leq 1}.
#' @param deviation factor of how much \eqn{n} can deviate from 0. Numeric value.
#' @param step size of the step for the search, between 0 and 1.
#' 
#' @return \code{sOpt} returns a numeric vector of the optimal exponent(s).
#' 
#' @examples 
#' sOpt(50, 0.6)
#' sOpt(100, 0.6)
#' 
#' @export
sOpt <- function(n, s = 0.6, deviation = 0, step = 0.001)
{
  if(deviation == 0) return(sOpt_c(n, s, step))
  
  sapply(n, function(n)
  {
    res <- sSizes(n, step = 0.01)
    # res$bn[which.min(abs(res$s[res$diff < 2] - s))]
    res$s[which.min(abs(res$s - 0.6) + deviation * res$diff)]
  })
}
