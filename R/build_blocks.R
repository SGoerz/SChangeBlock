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
#' for performing the \code{block_test}. 
#' 
#' @returns Strictly positive integer vector or matrix having the dimensions \code{n}. 
#' 
#' @export
build_blocks <- function(n, bn)
{
  ln <- floor(n / bn)
  bndiff <- n - ln * bn
  
  nblocks <- prod(bn)
  groups <- matrix(1:nblocks, nrow = bn[1])
  
  A <- kronecker(groups[0:bndiff[1], 0:bndiff[2]],
                 matrix(1, nrow = ln[1] + 1, ncol = ln[2] + 1))
  B <- kronecker(matrix(groups[0:bndiff[1], (bndiff[2] + 1):bn[2]], nrow = bndiff[1]),
                 matrix(1, nrow = ln[1] + 1, ncol = ln[2]))
  C <- kronecker(matrix(groups[(bndiff[1] + 1):bn[1], 0:bndiff[2]], nrow = bn[1] - bndiff[1]),
                 matrix(1, nrow = ln[1], ncol = ln[2] + 1))
  D <- kronecker(matrix(groups[(bndiff[1] + 1):bn[1], (bndiff[2] + 1):bn[2]], nrow = bn[1] - bndiff[1]),
                 matrix(1, nrow = ln[1], ncol = ln[2]))
  
  if(length(A) == 0) A <- NULL
  if(length(B) == 0) B <- NULL
  if(length(C) == 0) C <- NULL
  if(length(D) == 0) D <- NULL
  
  rbind(cbind(A, B), cbind(C, D))
}

sOpt2 <- function(n, s = 0.6)
{
  sapply(n, function(n)
  {
    res <- sSizes(n, step = 0.01)
    # res$bn[which.min(abs(res$s[res$diff < 2] - s))]
    res$bn[which.min(abs(res$s - 0.6) + 0.05 * res$diff)]
  })
}
