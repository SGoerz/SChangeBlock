test_that("outout from Mu has the correct format", {
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(length(Mu(X, c(2, 2))), 9)
  expect_equal(length(Mu(X, c(3, 3))), 4)
  # expect_equal(length(Mu(X, c(3, 2), c(1, 1))), 12)
})

test_that("Mu computes the correct block sums", {
  set.seed(1045)
  # vector:
  set.seed(1045)
  x <- rnorm(100)
  expect_equal(Mu(x, 50), c(mean(x[1:50]), mean(x[51:100])))
  
  # matrix:
  set.seed(1045)
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(Mu(X, c(2, 2)), c(0.6946552, -0.1936238, -0.4166834, 
                                 -0.4640347, -0.1491637, -0.112317, 
                                 0.08525641, -0.4745453, -0.5687769), 
               tolerance = 1e-7)
  expect_equal(Mu(X, c(3, 3)), c(-0.04461205, -0.327058, -0.4336385, 0.09453829), 
               tolerance = 1e-7)
  
  # errors:
  expect_error(Mu(1:3, 1))
  expect_error(Mu(rnorm(10), 0))
  expect_error(Mu(X, c(3, 0)))
})
