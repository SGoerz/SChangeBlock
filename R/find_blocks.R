#' Find blocks
#'
#' Find the blocks that contribute (most) to a high test statistic
#' 
#' @param x times series or random field to be investigated Either a numeric vector or a numeric matrix.
#' @param s parameter for the size of the blocks. 
#'          Either a single value 0 < s < 1 for determining the block length \eqn{l_n = [n^s]}, 
#'          or a vector with one entry per dimension between 0 and 1, 
#'          or a strictly positive integer vector or matrix indicating the group (or block) of the corresponding observation in X. 
#'          Default is \code{\link{sOpt}}\code{(n, 0.6)}.
#' @param fun Character string; one of "gmd" (default), "var", "varWithoutCentering".
#' @param u quantile after which blocks are marked significant.
#' 
#' @returns A numeric vector or matrix. Has the same dimensions and attributes as \code{x}.
#' 
#' @examples
#' # random field with a shift in the middle
#' X <- genField(c(50, 50), type = 2, H = 50)
#' Indices <- find_blocks(X, c(0.6, 0.6))
#' plot(Indices)
#' 
#' @export
find_blocks <- function(x, s, fun = "gmd", u = 0.95)
{
  attrX <- attributes(x)
  if(is(x, "vector") | is(x, "ts")) x <- matrix(x)
  n <- dim(x)
  
  if(missing(s)) s <- sOpt(n, 0.6)
  
  if(length(s) == 1 || length(s) == length(n))
  {
    if(any(s <= 0 | s >= 1)) stop("If s is meant as an exponent, it must be between 0 and 1.")
    ln <- round(n^s)
    lprod <- prod(ln)
    bn <- floor(n / ln)
    mu <- Mu(x, l = ln, includeEdge = TRUE)
  } else if( ((is(x, "vector") && is(s, "vector")) && length(x) == length(s)) || 
             (length(dim(x)) == length(dim(s)) && all(dim(x) == dim(s))) )
  {
    if(length(grep("var", fun)) < 1) stop("Individual block sizes only possible for 'var'")
    ln <- table(s)
    lprod <- ln
    bn <- length(ln)
    attributes(ln) <- NULL
    mu <- Mu(x, group = s)
  } else
  {
    stop("no valid option for s")
  }
  
  if(grepl("var", fun))
  {
    if(fun == "var")
    {
      S <- ((mu^2 - mean(mu)^2) * lprod / var(as.vector(x))) / sqrt(2)
    } else if(fun == "varWithoutCentering")
    {
      sigma <- sum(x^2) / (prod(n) - 1) 
      S <- ((mu^2 * lprod / sigma) - 1) / sqrt(2)
    }
    
    if(length(mu) == length(ln))
    {
      S <- matrix(S[s], ncol = ncol(s))
      index <- S > qnorm(u)
    } else
    {
      S <- matrix(S, ncol = bn[2])
      index <- as.logical(kronecker((S > qnorm(u)), matrix(1, nrow = ln[1], ncol = ln[2])))
    }
  } else if(fun == "gmd")
  {
    distMu <- dist(mu)
    S <- as.matrix(sqrt(prod(ln)) * distMu / var(as.vector(x)) - 2 / sqrt(pi))
    S <- S > qnorm(u, 0, 4 / 3 + 8 / pi * (sqrt(3) - 2))
    
    index <- integer(length(mu))
    
    for(i in 1:length(mu))
    {
      number <- max(colSums(S))
      if(number == 0) break
      
      temp <- which(colSums(S) == number)
      
      index[i] <- if(length(temp) > 1) sample(temp, 1) else temp
      
      S[index[i], ] <- 0
      S[, index[i]] <- 0
    }
    
    I <- matrix(0, nrow = bn[1], ncol = bn[2])
    I[index] <- 1
    index <- as.logical(kronecker(I, matrix(1, nrow = ln[1], ncol = ln[2])))
  }
  
  attributes(index) <- attrX
  return(index)
}
