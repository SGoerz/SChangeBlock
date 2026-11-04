test_that("bandwidth estimation returns the correct format", {
  x <- arima.sim(list(ar = 0.5), 100)
  X <- genField(c(50, 50), Theta = genTheta(10, 0.4))

  expect_equal(length(bandwidth(x, 0.3, 0.3)), 1)
  expect_equal(length(bandwidth(X, 0.3, 0.3)), 2)
})


test_that("de-correlation output has the correct format", {
  # random field
  X <- genField(c(20, 20), Theta = genTheta(10, 0.3))
  
  Y <- decorr(X, lags = bandwidth(X, 0.3, 0.3))
  expect_equal(class(X), class(Y))
  
  Y <- decorr(X, lags = bandwidth(X, 0.3, 0.3), separable = TRUE)
  expect_equal(class(X), class(Y))
  
  # time series
  x <- arima.sim(model = list(ar = 0.5), 100)
  
  y <- decorr(x, bandwidth(x, 0.3, 0.3))
  expect_equal(class(x), class(y))
})

test_that("de-correlation is working correctly", {
  # random field
  n <- c(20, 20)
  lags <- n^(1/3) * 0.9
  X <- genField(n, Theta = genTheta(10, 0.3))
  
  Y1 <- decorr(X, lags = n[1]^(1/3) * 0.9)
  Y2 <- decorr(X, lags = lags)
  
  expect_equal(Y1, Y2)
  
  ## centering not necessary if original random field is centered
  Y1 <- decorr(X, lags)
  Y2 <- decorr(X - mean(X), lags, withoutCentering = TRUE)
  
  expect_equal(Y1, Y2)
  
  ## separable autocovariance
  Y2 <- decorr(X, lags, separable = TRUE)
  
  # time series
  ## check if arima parameter is significant (to 10% level)
  x <- arima.sim(model = list(ar = 0.5), 100)
  y <- decorr(x, bandwidth(x, 0.3, 0.3))
  
  arima.x <- arima(x, c(1, 0, 0))
  arima.y <- arima(y, c(1, 0, 0))
  
  expect_lt((1 - pnorm(abs(arima.x$coef[1]) / sqrt(arima.x$var.coef[1])))*2, 0.1)
  expect_gt((1 - pnorm(abs(arima.y$coef[1]) / sqrt(arima.y$var.coef[1])))*2, 0.1)
})

