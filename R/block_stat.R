#' Block test statistic 
#'
#' Computes test statistics of the block test on structural changes.
#'
#' @param x times series or random field to be tested. Either a numeric vector or a numeric matrix.
#' @param s parameter for the size of the blocks. 
#'          Either a single value 0 < s < 1 for determining the block length \eqn{l_n = [n^s]}, 
#'          or a vector with one entry per dimension between 0 and 1, 
#'          or a strictly positive integer vector or matrix indicating the group (or block) of the corresponding observation in X. 
#'          Default is \code{\link{sOpt}}\code{(n, 0.6)}.
#' @param fun Character string; one of "gmd" (default), "var", "jb", "grubbs", "ANOVA".
#' @param varEstim variance estimator or variance estimation of the whole field or times series.
#'                 Either a function to estimate the variance with, or a numeric value.
#'                 
#' @details
#' First, the time series or random field is divided into blocks and the means of 
#' the blocks are computed. Then the function \code{fun} is applied to the block means: 
#' * gmd: Gini's mean difference
#' * var: Ordinary variance estimator
#' * jb: Jarque-Bera test
#' * grubbs: Grubbs test for outliers
#' * ANOVA: simple ANOVA.
#' 
#' @returns A numeric value. 
#' For \code{fun} = \code{"grubbs"} it has the attribute \code{n} 
#' indicating the number of blocks, i.e. the number of observations used in the Grubbs test.
#' For \code{fun} = \code{"ANOVA"} it has the attributes \code{k} (number of blocks) and
#' \code{N} (total number of observations).
#' 
#' @references \enc{Görz}{Goerz}, S. and Fried, R. (2025). "Detecting changes in the mean of spatial random fields on a regular grid", \emph{arXiv preprint} \href{https://arxiv.org/abs/2512.11599}{arXiv:2512.11599}
#' @references Schmidt, S. K. (2024). "Detecting changes in the trend function of heteroscedastic time series", Bernoulli, 30(4), 2598-2622.
#'  
#' @seealso [block_pValue]
#' 
#' @examples
#' # time series with a shift 
#' x <- arima.sim(model = list(ar = 0.5), n = 100)
#' x[1:50] <- x[1:50] + 1
#' block_stat(x, sOpt(100, 0.6))
#' 
#' # field without shift and ordinary variance
#' X <- genField(c(50, 50))
#' block_stat(X, sOpt(50, 0.6), "var")
#' 
#' # field with a shift and ordinary variance
#' X <- genField(c(50, 50), type = 2)
#' block_stat(X, sOpt(50, 0.6), "var")
#' 
#' # GMD test statistic, scaling variance estimated by the mad
#' block_stat(X, 0.6, fun = "var", varEstim = mad)
#'
#' @importFrom methods is
#' @export
block_stat <- function(x, s, fun = "gmd", varEstim = var)
{
  if(is.vector(x) | is.ts(x)) n <- length(x) else n <- dim(x)
  
  if(missing(s)) s <- sOpt(n, 0.6)

  if(length(s) == 1 || length(s) == length(n))
  {
    if(any(s <= 0 | s >= 1)) stop("If s is meant as an exponent, it must be between 0 and 1.")
    # CAREFUL: bn is the product of all dimensions, but ln is still a vector since
    # it is needed for the ANOVA
    ln <- round(n^s)
    bn <- prod(floor(n / ln))
    m <- Mu(x, l = ln)
  } else if( ((is(x, "vector") && is(s, "vector")) && length(x) == length(s)) || 
             (length(dim(x)) == length(dim(s)) && all(dim(x) == dim(s))) )
  {
    if(length(grep("var", fun)) < 1) stop("Individual block sizes only possible for 'var'")
    ln <- table(s)
    bn <- length(ln)
    attributes(ln) <- NULL
    m <- Mu(x, group = s)
  } else
  {
    stop("no valid option for s")
  }
  
  if(is.function(varEstim))
  {
    sigma_x <- varEstim(as.vector(x))
    sigma_mu <- varEstim(m)
  } else
  {
    sigma_x <- varEstim
    ### ????
    sigma_mu <- varEstim / prod(ln)
  }
  
  
  if(fun == "gmd")
  {
    # Gini's mean difference
    res <- sqrt(bn) * (sqrt(prod(ln)) / sqrt(sigma_x) * GMD(m) - 2 / sqrt(pi)) /
      sqrt(4/3 + 8/pi * (sqrt(3) - 2))
  } else if(fun == "var")
  {
    # variance
    # res <- (var(m * sqrt(prod(ln))) * (bn - 1) / sqrt(bn)) / (sigma_x * sqrt(2)) - (bn - 1) / sqrt(2 * bn)
    kappa <- 2#mean((m * sqrt(prod(ln)))^4) - 1
    if(length(ln) == length(m)) nuSq <- ln * m^2 else nuSq <- prod(ln) * m^2
    res <- (1 / sigma_x * (sum(nuSq) - prod(n) * mean(m)^2) - bn + 1) / sqrt(kappa * bn)
  } else if(fun == "varWithoutCentering")
  {
    sigma_x <- sum(as.vector(x)^2) / (prod(n) - 1)
    if(length(ln) == length(m)) nuSq <- ln * m^2 else nuSq <- prod(ln) * m^2
    res <- (sum(nuSq) - bn * sigma_x) / (sigma_x * sqrt(2 * bn))
  } else if(fun == "jb")
  {
    # Jarque-Bera
    s <- skewness(m)
    k <- kurtosis(m)
    
    res <- bn / 6 * (s^2 + (k - 3)^2 / 4)
  } else if(fun == "ks")
  {
    # KS statistics
    normvalues <- pnorm(sort(m), 0, sqrt(sigma_x / prod(ln)))
    res <- max(abs(c(normvalues - (1:bn) / bn, normvalues - (0:(bn-1)) / bn))) * sqrt(bn)
  } else if(fun == "grubbs")
  {
    # Grubbs test for outliers
    res <- grubbs(m, sigma_mu)
    attr(res, "n") <- bn
  } else if(fun == "ANOVA")
  {
    mu.bar <- mean(m)
    
    MST <- sigma_mu * prod(n)
    MSE <- sum(mu(x, ln, function(x) var(as.vector(x))))
    
    res <- MST / MSE
    #if(locfun == "median") res <- res * 2 / pi
    
    attr(res, "k") <- bn
    attr(res, "N") <- prod(n)
  } else
  {
    res <- m
  } 
  # else
  # {
  #   stop("Invalid argument for 'fun'")
  # }
  
  attr(res, "blocksize") <- ln
  
  return(res)
}
