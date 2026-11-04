test_that("cutoff output has the correct format", {
  # random field
  X <- genField(c(20, 20), Theta = genTheta(10, 0.3))
  
  Y <- cutoff(X, ln = c(4, 4), ex = c(1, 1))
  expect_equal(class(Y), class(X))
  expect_equal(dim(Y), dim(X) / 4)
  expect_equal(Y, cutoff(X[1:19, 1:19], ln = c(4, 4), ex = c(1, 1)))
  
  Y <- cutoff(X, ln = c(6, 6), ex = c(1, 1))
  expect_equal(dim(Y), floor(dim(X) / 6) + 1)
  
  # time series
  x <- arima.sim(model = list(ar = 0.5), 100)
  
  y <- cutoff(x, 4, ex = 1)
  expect_equal(class(y), class(x))
  expect_equal(tsp(y), c(1, (100 / 4 - 1) * 4 + 1, 1 / 4))
  
  x <- ts(x, start = 2, end = 200, frequency = 0.5)
  y <- cutoff(x, 4, ex = 1)
  expect_equal(tsp(y), c(2, (200 / 4 - 2) * 4 + 2, 0.5 / 4))
})

test_that("cutoff produces correct output", {
  # random field
  X <- genField(c(10, 10), Theta = genTheta(10, 0.3))
  Y <- cutoff(X, ln = c(4, 4), ex = c(1, 1))
  
  expect_equal(Y, matrix(c(mean(X[1:3, 1:3]), mean(X[5:7, 1:3]), mean(X[9:10, 1:3]), 
                           mean(X[1:3, 5:7]), mean(X[5:7, 5:7]), mean(X[9:10, 5:7]), 
                           mean(X[1:3, 9:10]), mean(X[5:7, 9:10]), mean(X[9:10, 9:10])), ncol = 3))
  
  # time series
  ## check if arima parameter is significant (to 10% level)
  set.seed(1)
  x <- arima.sim(model = list(ar = 0.5), 100)
  y <- cutoff(x, 3, 1)
  
  arima.x <- arima(x, c(1, 0, 0))
  arima.y <- arima(y, c(1, 0, 0))
  
  expect_lt((1 - pnorm(abs(arima.x$coef[1]) / sqrt(arima.x$var.coef[1])))*2, 0.1)
  expect_gt((1 - pnorm(abs(arima.y$coef[1]) / sqrt(arima.y$var.coef[1])))*2, 0.1)
})