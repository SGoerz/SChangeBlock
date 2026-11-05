test_that("skewness and kurtosis return the correct values", {
  expect_equal(skewness(rnorm(100000)), 0, tolerance = 0.05)
  expect_equal(kurtosis(rnorm(100000)), 3, tolerance = 0.05)
})

test_that("Grubbs test and corresponding critical values are correctly computed", {
  y <- replicate(10000, {x <- rnorm(100); grubbs(x)})
  expect_equal(mean(y > crit.grubbs(100, 0.05)), 0.05, tolerance = 0.05)
})

test_that("rmix returns the correct distribution", {
  set.seed(1895)
  x <- rmix(100, 0.1)
  set.seed(1895)
  y <- rnorm(100)
  index <- sample(100, rbinom(1, 100, 0.1))
  y[index] <- y[index] + 10
  
  expect_equal(x, y)
  
  set.seed(95)
  x <- rmix(100, q = 0.1, h = 0, sigma = 10)
  set.seed(95)
  y <- rnorm(100) 
  index <- sample(100, rbinom(1, 100, 0.1))
  y[index] <- y[index] * 10
  expect_equal(x, y)
})