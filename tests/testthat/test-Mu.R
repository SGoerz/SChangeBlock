test_that("output from Mu has the correct format", {
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(length(Mu(X, l = c(2, 2))), 9)
  expect_equal(length(Mu(X, l = c(3, 3))), 4)
  expect_equal(length(Mu(X, group = matrix(c(rep(1, 6), rep(2, 30)), ncol = 6))), 2)
  # expect_equal(length(Mu(X, c(3, 2), c(1, 1))), 12)
  
  x <- rnorm(6)
  expect_equal(length(Mu(x, l = 2)), 3)
})

test_that("Mu computes the correct block sums", {
  # vector:
  set.seed(1045)
  x <- rnorm(100)
  expect_equal(Mu(x, l = 50), c(mean(x[1:50]), mean(x[51:100])))
  expect_equal(Mu(x, group = c(rep(1, 10), rep(2, 90))), c(mean(x[1:10]), mean(x[11:100])))
  
  # matrix:
  set.seed(1045)
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(Mu(X, l = c(2, 2)), c(0.6946552, -0.4640347, 0.08525641,
                                     -0.1936238, -0.1491637, -0.4745453,
                                     -0.4166834, -0.112317, -0.5687769), 
               tolerance = 1e-7)
  expect_equal(Mu(X, l = c(3, 3)), c(-0.04461205, -0.4336385, -0.327058, 0.09453829), 
               tolerance = 1e-7)
  
  ## include edges:
  X <- matrix(rnorm(25), ncol = 5)
  
  expect_equal(Mu(X, l = c(1, 1)), as.vector(X))
  expect_equal(Mu(X, l = c(1, 1), includeEdge = TRUE), as.vector(X))
  
  X <- matrix(rnorm(30), ncol = 5)
  
  expect_equal(Mu(X, l = c(2, 2), includeEdge = FALSE), 
               c(mean(X[1:2, 1:2]), mean(X[3:4, 1:2]), mean(X[5:6, 1:2]), 
                 mean(X[1:2, 3:4]), mean(X[3:4, 3:4]), mean(X[5:6, 3:4])))
  
  expect_equal(Mu(X, l = c(2, 2), includeEdge = TRUE), 
               c(mean(X[1:2, 1:2]), mean(X[3:4, 1:2]), mean(X[5:6, 1:2]), 
                 mean(X[1:2, 3:4]), mean(X[3:4, 3:4]), mean(X[5:6, 3:4]), 
                 mean(X[1:2, 5]), mean(X[3:4, 5]), mean(X[5:6, 5])))
  
  X <- matrix(rnorm(30), ncol = 6)
  
  expect_equal(Mu(X, l = c(2, 2), includeEdge = FALSE), 
               c(mean(X[1:2, 1:2]), mean(X[3:4, 1:2]), 
                 mean(X[1:2, 3:4]), mean(X[3:4, 3:4]), 
                 mean(X[1:2, 5:6]), mean(X[3:4, 5:6])))
  
  expect_equal(Mu(X, l = c(2, 2), includeEdge = TRUE), 
               c(mean(X[1:2, 1:2]), mean(X[3:4, 1:2]), mean(X[5, 1:2]), 
                 mean(X[1:2, 3:4]), mean(X[3:4, 3:4]), mean(X[5, 3:4]), 
                 mean(X[1:2, 5:6]), mean(X[3:4, 5:6]), mean(X[5, 5:6])))
  
  # errors:
  expect_error(Mu(1:3, group = 1))
  expect_error(Mu(rnorm(10), l = 0))
  expect_error(Mu(X, l = c(3, 0)))
  expect_error(Mu(c(TRUE, FALSE), l = 1))
  expect_error(Mu(1:2))
})
