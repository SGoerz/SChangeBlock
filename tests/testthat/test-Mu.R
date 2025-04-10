test_that("outout from Mu has the correct format", {
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(length(Mu(X, c(2, 2))), 9)
  expect_equal(length(Mu(X, c(3, 3))), 4)
  # expect_equal(length(Mu(X, c(3, 2), c(1, 1))), 12)
})

test_that("Mu computes the correct block sums", {
  set.seed(1045)
  X <- matrix(rnorm(36), ncol = 6)
  
  expect_equal(Mu(X, c(2, 2)), c(0.6946552, -0.1936238, -0.4166834, 
                                 -0.4640347, -0.1491637, -0.112317, 
                                 0.08525641, -0.4745453, -0.5687769), 
               tolerance = 1e-7)
  expect_equal(Mu(X, c(3, 3)), c(-0.04461205, -0.327058, -0.4336385, 0.09453829), 
               tolerance = 1e-7)
  # expect_equal(Mu(X, c(3, 2), c(1, 1)), c(1.3644447, 0.05410421, 0.09391388, 
  #                                         -0.2857292, -0.4905827, -0.6359239, 
  #                                         -0.08935994, 0.1169215, 0.06596551, 
  #                                         -0.32818816, -0.2968408, -0.830319))
})
