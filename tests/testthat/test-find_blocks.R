test_that("find_blocks output has the correct format", {
  X <- genField(c(20, 20))
  Y <- find_blocks(X, sOpt(c(20, 20)), fun = "var", u = 0.8)
  
  expect_equal(attributes(X), attributes(Y))
  
  x <- arima.sim(list(ma = 0.2), 100)
  y <- find_blocks(x, 0.646)
  
  expect_equal(attributes(x), attributes(y))
  
  expect_error(find_blocks(X, s = 2))
  expect_error(find_blocks(X, s = c("a", "b", "c")))
  expect_error(find_blocks(X, s = build_blocks(c(20, 20), s = 0.6), fun = "gmd"))
})

test_that("find_blocks produces the correct output", {
  n <- c(50, 50)
  X <- genField(n, type = 2, H = 50)
  
  ln <- floor(n^0.6)
  bn <- round(n / ln)
  mu <- Mu(X, l = ln)
  
  # var 
  Y <- find_blocks(X, 0.6, fun = "var", u = 0.8)
  
  S <- prod(ln) * (mu^2 - mean(X)^2) / (sqrt(2) * var(as.vector(X))) > qnorm(0.8)
  S <- as.logical(kronecker(matrix(S, nrow = bn[1]), matrix(1, nrow = ln[1], ncol = ln[2])))
  attributes(S) <- attributes(X)
  
  expect_equal(S, Y)
  
  # var without centering
  s <- build_blocks(n, s = 0.6)
  Y <- find_blocks(X, s, fun = "varWithoutCentering", u = 0.9)
  
  sigma <- sum(X^2) / (prod(n) - 1)
  S <- (prod(ln) * mu^2 / sigma - 1) / sqrt(2) > qnorm(0.9)
  S <- as.logical(kronecker(matrix(S, nrow = bn[1]), matrix(1, nrow = ln[1], ncol = ln[2])))
  attributes(S) <- attributes(X)
  
  expect_equal(S, Y)
  
   # gmd
  set.seed(110)
  Y <- find_blocks(X, fun = "gmd", u = 0.95)
  
  set.seed(110)
  S <- as.matrix(dist(mu) * sqrt(prod(ln)) / var(as.vector(X)) - 2 / sqrt(pi))
  S <- S > qnorm(0.95, 0, 4 / 3 + 8 / pi * (sqrt(3) - 2))
  
  index <- numeric(ncol(S))
  
  for(i in seq_along(index))
  {
    candidate <- max(colSums(S))
    if(candidate == 0) break
    candidate <- which(colSums(S) == candidate)
    
    index[i] <- if(length(candidate) == 1) candidate else sample(candidate, 1)
    S[index[i], ] <- FALSE
    S[, index[i]] <- FALSE
  }
  
  S <- logical(length(index))
  S[index] <- TRUE
  
  S <- as.logical(kronecker(matrix(S, nrow = bn[1]), matrix(1, nrow = ln[1], ncol = ln[2])))
  attributes(S) <- attributes(X)
  
  expect_equal(S, Y)
})
