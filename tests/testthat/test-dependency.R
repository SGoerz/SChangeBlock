test_that("dependency matrix Phi is correct", {
  phi <- genPhi(1, 0.4)
  
  expect_true("matrix" %in% class(phi))
  expect_equal(dim(phi), c(3, 3))
  expect_equal(phi, matrix(c(0.4^2, 0.4, 0.4^2, 
                             0.4, 1, 0.4, 
                             0.4^2, 0.4, 0.4^2), ncol = 3))
  
  phi <- genPhi(2, c(0.8, 0.6, 0.4, 0.2))
  expect_equal(dim(phi), c(5, 5))
  expect_equal(phi, matrix(c(0.2, 0.4, 0.6, 0.4, 0.2, 
                             0.4, 0.6, 0.8, 0.6, 0.4, 
                             0.6, 0.8, 1, 0.8, 0.6, 
                             0.4, 0.6, 0.8, 0.6, 0.4, 
                             0.2, 0.4, 0.6, 0.4, 0.2), ncol = 5))
})

test_that("dependency is correctly added to the data", 
{
  E <- matrix(rnorm(25), ncol = 5)
  phi <- genPhi(2, 0.4)
  X1 <- dependency(E, q_ = 2, param_ = 0.4)
  X2 <- dependency(E, Phi_ = phi)
  Xtest <- sum(E * phi)
  
  expect_equal(X1, X2)
  expect_true("matrix" %in% class(X1))
  expect_equal(dim(X1), c(1, 1))
  expect_equal(X1[1], Xtest)
  expect_equal(X2[1], Xtest)
})